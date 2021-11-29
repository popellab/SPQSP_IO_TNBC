#ifndef __OUTPUT_HUB__
#define __OUTPUT_HUB__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using std::ofstream;
using std::string;

//! Handles all output stream to files
class FileOutputHub
{
public:
	FileOutputHub(void);
	~FileOutputHub();

	//! initialize output overseer. Need to be called in main before using the object.
	void setup(long seed, string baseDirectory);

	//! get pointer to log file streaml
	ofstream& getLogFstream() { return _log; };
	//! get pointer to general stats file stream.
	ofstream& getStatsFstream(bool, int);
	//! get reference to QSP output file stream;
	ofstream& get_lymph_blood_QSP_stream(void) { return _lymph_blood_QSP; };

	//! create a new output stream for grid snap shot and return the pointer
	ofstream& getNewGridToSnapshotStream(unsigned long time, string tag );
	//! create a new output stream for serilization and return the pointer
	ofstream& getNewSaveStateStream(unsigned long time, string tag );
	//! create a new output stream for other purposes.
	ofstream& getNewMakeshiftStream(string tag);
	
private:

	//! seed used for simulation
	long _seed;
	//! base directory path for all output
	string _outDirBase; //!< end without "/"
	//! log output stream
	ofstream _log;
	//! general stats output stream
	std::vector<ofstream> _generalStats_core;
	std::vector<ofstream> _generalStats_margin;
	//! Lymph-Blood QSP ode
	ofstream _lymph_blood_QSP;
	//! grid snapshot output stream
	ofstream _gridSnapshotFstream;
	//! serialization output stream
	ofstream _saveStateStream;
	//! makeshift file stream
	ofstream _makeshiftSteam;

};

#endif
