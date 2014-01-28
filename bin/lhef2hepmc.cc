// -*- C++ -*-

// lhef2hepmc: a Les Houches Event Format (LHEF) -> HepMC IO_GenEvent
// MC generator event file converter
// Author: Andy Buckley, Les Houches 2011
// Hat tip: Leif Lonnblad for writing the LHEF parser that actually makes it possible!

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "Bianchi/TTHStudies/interface/LHEF.h"
#include <cstdlib>
#include <string>
#include <cassert>
#include <iostream>

#include "TSystem.h"
#include "TROOT.h"

using namespace std;
using namespace HepMC;

// In case IDWTUP=+/-4 one has to keep track of the accumulated weights and event numbers
// to evaluate the cross section on-the-fly. The last evaluation is the one used.
// Better to be sure that crossSection() is never used to fill the histograms, but only in the finalization stage, by reweighting the histograms with crossSection()/sumOfWeights()

double accumulated_weight = 0.0;
double accumulated_weight_squared = 0.0;
int event_number = 0;


int main(int argc, char** argv) {

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  std::cout << "lhe --> HepMC converter" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  // Look for a help argument
  for (int i = 1; i < argc; ++i) {
    const string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      cout << argv[0] << ": an LHEF to HepMC event format converter" << endl;
      cout << "Usage: " << argv[0] << "[<infile> <outfile>]" << endl;
      cout << "If called with no args, reads from stdin and writes to stdout." << endl;
      cout << "Use filename '-' to explicitly specify reading/writing on stdin/out." << endl;
      exit(0);
    }
  }


  // Choose input and output files from the command line
  string infile("-"), outfile("-");
  if (argc == 3) {
    infile = argv[1];
    outfile = argv[2];
  } else if (argc == 2 || argc > 3) {
    cerr << "Usage: " << argv[0] << "[<infile> <outfile>]" << endl;
    exit(1);
  }


  // Init readers and writers
  LHEF::Reader* reader = 0;
  if (infile == "-") {
    reader = new LHEF::Reader(cin);
  } else {
    reader = new LHEF::Reader(infile);
  }
  IO_GenEvent* writer = 0;
  if (outfile == "-") {
    writer = new IO_GenEvent(cout);
  } else {
    writer = new IO_GenEvent(outfile);
  }

  if( reader==0 ) cout << "No reader" << endl;

  // Event loop
  while (reader->readEvent()) {
    GenEvent evt;
    evt.use_units(Units::GEV, Units::MM);
    const double weight = reader->hepeup.weight();
    evt.weights().push_back(weight);

    // Cross-section
    #ifdef HEPMC_HAS_CROSS_SECTION
    HepMC::GenCrossSection xsec;
    const int idwtup = reader->heprup.IDWTUP;

    double xsecval = -1.0;
    double xsecerr = -1.0;
    if (abs(idwtup) == 3) {
      xsecval = reader->heprup.XSECUP[0];
      xsecerr = reader->heprup.XSECUP[1];
      //cout << "Read cross-section = " << xsecval << " +- " << xsecerr << endl;
    } else if (abs(idwtup) == 4) {
      accumulated_weight += weight;
      accumulated_weight_squared += weight*weight;
      event_number += 1;
      xsecval = accumulated_weight/event_number;
      // xsecerr = sqrt((accumulated_weight_squared/event_number - xsecval*xsecval)/event_number);
      double xsecerr2 = (accumulated_weight_squared/event_number - xsecval*xsecval)/event_number;
      if (xsecerr2 < 0) {
        cerr << "WARNING: xsecerr^2 < 0, forcing to zero : " << xsecerr2 << endl;
        xsecerr2 = 0.0;
      }
      xsecerr = sqrt(xsecerr2);
      // cout << "Estimated cross-section = " << xsecval << " +- " << xsecerr << endl;
    } else  {
      cerr << "IDWTUP = " << idwtup << " value not handled yet. Stopping " << endl;
      exit(-1);
    }
    xsec.set_cross_section(xsecval, xsecerr);
    evt.set_cross_section(xsec);
    #endif

    GenVertex* v = new GenVertex();
    evt.add_vertex(v);
    FourVector beam1(0, 0, reader->heprup.EBMUP.first, reader->heprup.EBMUP.first);
    GenParticle* gp1 = new GenParticle(beam1, reader->heprup.IDBMUP.first, 4);
    v->add_particle_in(gp1);
    FourVector beam2(0, 0, -reader->heprup.EBMUP.second, reader->heprup.EBMUP.second);
    GenParticle* gp2 = new GenParticle(beam2, reader->heprup.IDBMUP.second, 4);
    v->add_particle_in(gp2);
    evt.set_beam_particles(gp1, gp2);

    for (int i = 0; i < reader->hepeup.NUP; ++i) {
      if (reader->hepeup.ISTUP[i] != 1) {
        //cout << reader->hepeup.ISTUP[i] << ", " << reader->hepeup.IDUP[i] << endl;
        continue;
      }
      FourVector p(reader->hepeup.PUP[i][0], reader->hepeup.PUP[i][1],
                   reader->hepeup.PUP[i][2], reader->hepeup.PUP[i][3]);
      GenParticle* gp = new GenParticle(p, reader->hepeup.IDUP[i], 1);
      gp->set_generated_mass(reader->hepeup.PUP[i][4]);
      v->add_particle_out(gp);
    }

    //evt.print();
    writer->write_event(&evt);
  }

  delete reader;
  delete writer;

  return EXIT_SUCCESS;
}
