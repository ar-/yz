/*
 * main.cxx
 * 
 * Copyright 2014 Andreas Altair Redmer <altair.ibn.la.ahad.sy@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include "main.h"
#include <algorithm>

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

void printHelp()
{
	cout << "help";
}

int main(int argc, char **argv)
{

	return GeneralPacker::pg_main(argc, argv);

	if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        cout << "help" << endl;
    }

    char * filename = getCmdOption(argv, argv + argc, "-f");

    if (filename)
    {
        // Do interesting things
        // ...
    }
    
	if(cmdOptionExists(argv, argv+argc, "--cpu"))
    {
        cout << "cpu" << endl;
    }
    

	cout << "dd" << endl;
	return 0;
}

