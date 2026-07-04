// Setup.h is part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// It is modified from the PYTHIA event generator, which is
// Copyright (C) 2008 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the settings database; and for error statistics.
// Flag: helper class with bool flags.
// Mode: helper class with int modes.
// Parm: (short for parameter) helper class with double parameters.
// Word: helper class with string words.
// Settings: maps of flags, modes, parms and words with input/output.

#ifndef Setup_H
#define Setup_H

#include <iostream>
#include "Information.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

class Setup {



public:

  // Constructor.
  Setup() {}

  // Initialize static pointer.
  static void initPtr(Info* infoPtrIn) {infoPtr = infoPtrIn;}
 
  // Read in database from specific file.
  static bool init(string startFile = "Parameters.xml", 
    bool append = false, ostream& os = cout) ;

  // Read in one update from a single line.
  static bool readString(string line, bool warn = true, 
    ostream& os = cout) ; 
 
  // Print out table of database, either all or only changed ones,
  // or ones containing a given string.
  static void listAll(ostream& os = cout) { 
    list( true, false, " ", os); } 
  static void listChanged(ostream& os = cout) { 
    list (false, false, " ", os); } 
  static void list(string match, ostream& os = cout) { 
    list (false, true, match, os); } 

  // Query existence of an entry.
  static bool isFlag(string keyIn) {
    return (flags.find(toLower(keyIn)) != flags.end()); }
  static bool isMode(string keyIn) { 
    return (modes.find(toLower(keyIn)) != modes.end()); }
  static bool isParm(string keyIn) {
    return (parms.find(toLower(keyIn)) != parms.end()); }
  static bool isWord(string keyIn) {
    return (words.find(toLower(keyIn)) != words.end()); }

  // Add new entry.
  static void addFlag(string keyIn, bool defaultIn) {
    flags[toLower(keyIn)] = Flag(keyIn, defaultIn); }  
  static void addMode(string keyIn, int defaultIn, bool hasMinIn, 
    bool hasMaxIn, int minIn, int maxIn) { modes[toLower(keyIn)] 
    = Mode(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }      
  static void addParm(string keyIn, double defaultIn, bool hasMinIn, 
    bool hasMaxIn, double minIn, double maxIn) { parms[toLower(keyIn)] 
    = Parm(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }  
  static void addWord(string keyIn, string defaultIn) {
    words[toLower(keyIn)] = Word(keyIn, defaultIn); }  
 
  // Give back current value, with check that key exists. 
  static bool   flag(string keyIn);
  static int    mode(string keyIn);
  static double parm(string keyIn);
  static string word(string keyIn); 
  
  // Change current value, respecting limits.
  static void flag(string keyIn, bool nowIn); 
  static void mode(string keyIn, int nowIn);
  static void parm(string keyIn, double nowIn); 
  static void word(string keyIn, string nowIn); 

  // Change current value, disregarding limits.
  static void forceMode(string keyIn, int nowIn) { 
    if (isMode(keyIn)) modes[toLower(keyIn)].valNow = nowIn; }
  static void forceParm(string keyIn, double nowIn) { 
    if (isParm(keyIn)) parms[toLower(keyIn)].valNow = nowIn; }
     
private:

  // Pointer to various information on the generation.
  static Info* infoPtr;

  // Map for bool flags.
  static map<string, Flag> flags;

  // Map for integer modes.
  static map<string, Mode> modes;

  // Map for double parms.
  static map<string, Parm> parms;

  // Map for string words.
  static map<string, Word> words;

  // Flag that initialization has been performed.
  static bool isInit;

  // Print out table of database, called from listAll and listChanged.
  static void list(bool listAll, bool listString, string match,
    ostream& os = cout) ; 

  // Useful functions for string handling.
  static string toLower(const string& name);
  static bool   boolString(string tag);
  static string attributeValue(string line, string attribute);
  static bool   boolAttributeValue(string line, string attribute);
  static int    intAttributeValue(string line, string attribute);
  static double doubleAttributeValue(string line, string attribute);

};

//**************************************************************************


#endif // Setup_H
