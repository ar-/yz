/*
 * packer.h
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

#ifndef packer_h
#define packer_h

#define kvp string,string
#define attrs map<kvp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cerrno>
#include <map>
#include <vector>

using namespace std;


namespace yztools
{
	class Packer
	{
	private:
		std::string getFileContentsString(const char *filename);
	public:
		virtual ~Packer() {}
		virtual void getAttribute(string &s) = 0;
		virtual attrs* getAttributes(string &s) = 0;
		
	};
}

#endif
