#ifndef FILEREADER_HH
#define FILEREADER_HH
#include "Particle.hpp"
#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <vector>


class FileReader
{

    public:

	FileReader()
	{
	}

	//register a new parameter with name key and initial int value
	void RegisterIntParameter(const std::string& key, int init)
	{
		paramTypes.insert(std::pair<std::string,std::string>(key,"INT"));
		intParams.insert(std::pair<std::string,int>(key,init));
	}

	//register a new parameter with name key and initial double value
	void RegisterDoubleParameter(const std::string& key, double init)
	{
		paramTypes.insert(std::pair<std::string,std::string>(key,"DOUBLE"));
		doubleParams.insert(std::pair<std::string,double>(key,init));
	}

	//register a new parameter with name key and initial string value
	void RegisterStringParameter(const std::string& key, const std::string &init)
	{
		paramTypes.insert(std::pair<std::string,std::string>(key,"STRING"));
		stringParams.insert(std::pair<std::string,std::string>(key,init));
	}

	//set a value for the key string with value in
	void SetParameter(const std::string &key, const std::string &in)
	{
		std::map<std::string,std::string>::iterator it=stringParams.find(key);
		it->second=in;
	}

	//set a value for the key string with value in
	void SetParameter(const std::string &key, double in)
	{
		std::map<std::string,double>::iterator it=doubleParams.find(key);
		it->second=in;
	}
	//set a value for the key string with value in
	void SetParameter(const std::string &key, int in)
	{
		std::map<std::string,int>::iterator it=intParams.find(key);
		it->second=in;
	}

	// get the int value of key 
	inline int GetIntParameter(const std::string &key) const
	{
		std::map<std::string,int>::const_iterator it=intParams.find(key);
		return it->second;
	}

	// get the double value of key 
	inline double GetDoubleParameter(const std::string &key) const
	{
		std::map<std::string,double>::const_iterator it=doubleParams.find(key);
		return it->second;
	}

	// get the string value of key 
	inline std::string GetStringParameter(const std::string &key) const
	{
		std::map<std::string,std::string>::const_iterator it=stringParams.find(key);
		return it->second;
	}
	
	//Tells which type is it
	inline std::string GetType(const std::string &key) const
	{
		std::map<std::string,std::string>::const_iterator it=paramTypes.find(key);
		if(it==paramTypes.end())
		{
			std::cout<<"ERROR: such "<<key<< " key was not found in the registerd keys\n";
			//return " ";
		}
		return it->second;
	
		
	}
	
	//try to read all registered parameters from file name
	bool ReadFile(const std::string &name)
	{
		std::ifstream input( name.c_str() );
		//std::cout<<input;
		std::string line;
		int state=0;
		if (input!=NULL) state=1;//Check if the file is good
		while( getline( input, line))
		{	if (!line.empty())
			{
				std::istringstream iss(line);
				std::string sub1,sub2;
        		iss >>std::ws>>sub1>>std::ws>>sub2;
    			if(sub1.length()!=0)
				{
    				std::string type=GetType(sub1);
    				if(type=="INT"){
    					int value=atoi(sub2.c_str());
    					//std::cout<<sub1<<"   "<<value<<"\n";
						SetParameter(sub1, value);
    				}
    				if(type=="DOUBLE"){
    					double valued=atof(sub2.c_str());
    					//std::cout<<sub1<<"   "<<value<<"\n";
						SetParameter(sub1, valued);
    				}
    				if(type=="STRING"){
    					//std::cout<<sub1<<"   "<<sub2<<"\n";
						SetParameter(sub1, sub2);
    				}
    					
    			}
    		}
		}
		if (state==1 )return true;
		  else return false;	
	}
		
   std::vector<Particle> pvec;
	
   //try to read DATA file
   bool read_data(const std::string &name)
   {
        int state=0;
	std::ifstream input( name.c_str() );
	//std::cout<<input;
	std::string line;
	
	if (input!=NULL) state=1;//Check if the file is good
	while( getline( input, line))
	{	if (!line.empty())
		{
			std::istringstream iss(line);
			double  sub1,sub2,sub3,sub4,sub5,sub6,sub7;
        		iss >>sub1>>sub2>>sub3>>sub4>>sub5>>sub6>>sub7;
			//pvec.push_back(Particle(atof(sub1.c_str()), atof(sub2.c_str()), atof(sub3.c_str()), atof(sub4.c_str()), atof(sub5.c_str()), atof(sub6.c_str()), atof(sub7.c_str())));
			pvec.push_back(Particle(-1,sub1,sub2,sub3,sub4,sub5,sub6,sub7,0,0,0));
    		}
	}
	
	if (state==1 )return true;
	else return false;
	
    }

	//print out all parameters to std:out
	void PrintParameters() const
	{
		std::map<std::string,int>::const_iterator iti=intParams.begin(), end_iti=intParams.end();
		std::map<std::string,double>::const_iterator itd=doubleParams.begin(), end_itd=doubleParams.end();
		std::map<std::string,std::string>::const_iterator its=stringParams.begin(), end_its=stringParams.end();
		for(;iti!=end_iti;++iti)
			std::cout<<iti->first<<"  "<<iti->second<<"\n";
		for(;itd!=end_itd;++itd)
			std::cout<<itd->first<<"  "<<itd->second<<"\n";
		for(;its!=end_its;++its)
			std::cout<<its->first<<"  "<<its->second<<"\n";
	}

    private:
	std::map<std::string,std::string> paramTypes;
	std::map<std::string,int> intParams;
	std::map<std::string,double> doubleParams;
	std::map<std::string,std::string> stringParams;
    
};

#endif //FILEREADER_HH

