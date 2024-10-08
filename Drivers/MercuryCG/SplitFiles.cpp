//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

///takes data and fstat files and splits them into *.data.???? and *.fstat.???? files

#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h> 
#include <cstdio>
#include <cstdlib>
#include <Logger.h>


class CFile {

public:

	///Constructor
	explicit CFile(std::string name) {
		//set file names
		data_filename.str("");
		data_filename << name << ".data";
		fstat_filename.str("");
		fstat_filename << name << ".fstat";
		
		//open in-streams
		data_file.open(data_filename.str().c_str(), std::fstream::in);
		fstat_file.open(fstat_filename.str().c_str(), std::fstream::in);

		if (data_file.fail())
		{
            logger(ERROR, "Input file % not found", data_filename.str());
		} else
        {
            if (fstat_file.fail())
            {
                logger(WARN, "WARN: Input file % not found; only data files will be processed", fstat_filename.str()
                );
                logger(INFO, "File opened: %", data_filename.str());
            }
            else
            {
                logger(INFO, "Files opened: % and %", data_filename.str(), fstat_filename.str());
            }
        }

    }

	///Destructor
	~CFile()
    {
        data_file.close();
        fstat_file.close();
        logger(INFO, "Files closed: % and %", data_filename.str(), fstat_filename.str());
    }
		
	bool copy(unsigned int stepsize, unsigned int counter) {
		//return copy_last_time_step();
        if (fstat_file.fail())
            return copy_data(stepsize,counter);
        else
            return copy_data(stepsize,counter) && copy_fstat(stepsize,counter);
	}

	bool copy_data(unsigned int stepsize, unsigned int counter) {
		unsigned int N;
		std::string line;
		std::stringstream output_filename;
		std::fstream output_file;				
		
		data_file >> N;
		while (data_file.good()) {
			//set output_filename
			std::stringstream output_filename("");
			output_filename << data_filename.str() << ".";
            if (counter < 1000) output_filename << "0";
            if (counter < 100) output_filename << "0";
            if (counter < 10) output_filename << "0";
            output_filename << counter;
            counter++;
            //logger(INFO, "Outputfile: %", output_filename.str());
            
            //open, write, close output file
            output_file.open(output_filename.str().c_str(), std::fstream::out);
            getline(data_file, line);
            output_file << N << line << std::endl;
            logger(INFO, "%%", N, line);
            for (unsigned int i = 0; i < N; i++)
            {
                getline(data_file, line);
                output_file << line << std::endl;
            }
            output_file.close();
            double doubleN = -1;
            data_file >> doubleN;
            N = doubleN;
            
            // converting N to an integer; skipping the line if there is a problem (this happens when there is a corrupt data file)
			while (data_file.good() && doubleN != N)
            {
                getline(data_file, line, '\n');
                logger(WARN, "Skipping bad line in data file: %", doubleN);
                data_file >> doubleN;
                N = doubleN;
            }

            //step over some time steps
			for(unsigned int j=1; j<stepsize; j++) {
				for (unsigned int i=0; i<N+1; i++) getline(data_file,line);
				data_file >> N;
			}

		}
		return true;
	}

	bool copy_fstat(unsigned int stepsize, unsigned int counter) {
		std::string line;
		std::stringstream output_filename;
		std::fstream output_file;				
		
		getline(fstat_file,line);
		while (fstat_file.good()) {
            //set output_filename
            output_filename.str("");
            output_filename << fstat_filename.str() << ".";
            if (counter < 1000) output_filename << "0";
            if (counter < 100) output_filename << "0";
            if (counter < 10) output_filename << "0";
            output_filename << counter;
            counter++;
            
            //open, write, close output file
            output_file.open(output_filename.str().c_str(), std::fstream::out);
            logger(INFO, "%", line);
            output_file << line << std::endl;
            getline(fstat_file, line);
            output_file << line << std::endl;
            getline(fstat_file, line);
            output_file << line << std::endl;
            getline(fstat_file, line);
            while (line.c_str()[0] != '#' && fstat_file.good())
            {
                getline(fstat_file, line);
                output_file << line << std::endl;
            }
            output_file.close();

			//step over some time steps
			for(unsigned int j=1; j<stepsize; j++) {
				getline(fstat_file,line); 
				getline(fstat_file,line); 
				getline(fstat_file,line);
				while (line.c_str()[0] != '#'&&fstat_file.good()) {
					getline(fstat_file,line);
				}
			}
		}
		return true;
	}


private:
	///These store the save file names, 
	std::stringstream data_filename;
	std::stringstream fstat_filename;

	///Stream used for data files
	std::fstream data_file;
	std::fstream fstat_file;
};

int main(int argc, char *argv[])
{
	if (argc<2) {
        logger(WARN, "split_files problem_name [stepsize [initial_counter]]");
		return -1;
	}
	std::string name(argv[1]);
    logger(INFO, "Name: %", name);
    
    unsigned int stepsize = 1;
    if (argc > 2) stepsize = static_cast<unsigned int>(atoi(argv[2]));
    
    //defines the initial counter
    unsigned int counter = 0;
    if (argc > 3) counter = static_cast<unsigned int>(atoi(argv[3]));
    
    CFile files(name);
    files.copy(stepsize, counter);
    logger(INFO, "finished writing split files: %", name);
    return 0;
}
