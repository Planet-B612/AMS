#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif
#include "serialize.h"
#include <memory>
#include <fstream>

using namespace std;
class FileCtrl
{
public:
	static void mkdir_absence(const char* outFolder)
	{
#if defined(_WIN32)
		CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
		mkdir(outFolder, 0733); // can be used on non-Windows
#endif
	}

	/// Save a serialized file
	template <class T>
	static void save_file(const string filename, const T& output)
	{
		ofstream outfile(filename, std::ios::binary);
		if (!outfile.eof() && !outfile.fail())
		{
			StreamType res;
			serialize(output, res);
			outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
			outfile.close();
			res.clear();
			cout << "Save file successfully: " << filename << '\n';
		}
		else
		{
			cout << "Save file failed: " + filename << '\n';
			exit(1);
		}
	}

	/// Load a serialized file
	template <class T>
	static void load_file(const string filename, T& input)
	{
		ifstream infile(filename, std::ios::binary);
		if (!infile.eof() && !infile.fail())
		{
			infile.seekg(0, std::ios_base::end);
			const std::streampos fileSize = infile.tellg();
			infile.seekg(0, std::ios_base::beg);
			vector<uint8_t> res(fileSize);
			infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
			//infile.clear(); // added by me
			infile.close();
			input.clear();
			auto it = res.cbegin();
			//cout<<res.cend()-res.cbegin()<<", "<<sizeof(T)<<endl;
			input = deserialize<T>(it, res.cend());
			res.clear();
		}
		else
		{
			cout << "Cannot open file: " + filename << '\n';
			exit(1);
		}
	}

	/// Save graph structure to a file
	static void save_graph_struct(const string graphName, const Graph& vecGraph, const bool isReverse)
	{
		string postfix = ".vec.graph";
		if (isReverse) postfix = ".vec.rvs.graph";
		const string filename = graphName + postfix;
		save_file(filename, vecGraph);
	}

	/// Load graph structure from a file
	static void load_graph_struct(const string graphName, Graph& vecGraph, const bool isReverse)
	{
		string postfix = ".vec.graph";
		if (isReverse) postfix = ".vec.rvs.graph";
		const string filename = graphName + postfix;
		load_file(filename, vecGraph);
	}

	/// Get out-file name
	static string get_out_file_name(const string graphName, const string algName, const double Q)
	{
		return graphName + "_" + algName + "_Q" + to_string(Q);
	}

};

using TIO = FileCtrl;
using PIO = std::shared_ptr<FileCtrl>;