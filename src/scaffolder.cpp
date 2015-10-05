#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <stack>
#include <limits>
#include <boost/lexical_cast.hpp>
#include <fcntl.h>

using namespace std;
using namespace boost;

#ifndef ID_TYPE
#error define ID_TYPE 1 .. string or 2 .. int
#endif

#if ID_TYPE == 1
	#define ID_T string
	#define ID_REF_T string&
#else
	#define ID_T int
	#define ID_REF_T int
#endif

void get_contig_sizes(const string &infile,map<string,size_t> &result);

void usage(char *s)
{	cout << s << ": Usage:" << endl;	
	cout << "	-h		this help" << endl;
	cout <<	"	-s graphfile" << endl;
	cout << "	-f fastafile" << endl;
	cout << "	-o outfile" << endl;
	cout << "	-v		increase verbosity (multiple possible)" << endl;
	cout << endl;
}

class edge
{	public:
	ID_T	n1,n2;
	edge()
	:n1(),n2()
	{
	}
	edge(const ID_REF_T s1,const ID_REF_T s2)
	{	if (s1<s2)
		{	n1=s1;
			n2=s2;
		}
		else
		{	n2=s1;
			n1=s2;
		}	
		if (s1==s2)
		{	throw runtime_error("self edge ???");
		}
	}
	bool operator < (const edge &e) const
	{	if (n1!=e.n1)
		{	return n1<e.n1;
		}
		return n2<e.n2;
	}
	bool operator == (const edge &e)
	{	return n1==e.n1 && n2==e.n2;
	}
};	

class edge_data
{	public:
	int count;
	edge_data()
	{	count=0;
	}
	edge_data(int c)
	{	count=c;
	}
};	

ostream& operator << (ostream &o,const edge &e)
{	return o << e.n1 << '-' << e.n2;
}

int get_count_edges(const map<ID_T,list<ID_T> > &c)
{	int ret=0;
	for (map<ID_T,list<ID_T> >::const_iterator i=c.begin(); i!=c.end(); ++i)
	{	ret+=i->second.size();
	}
	return ret/2;
}	

void remove_from_connections(map<ID_T,list<ID_T> > &connections,const edge &e)
{	list<ID_T>::iterator i=find(connections[e.n1].begin(),connections[e.n1].end(),e.n2);
	connections[e.n1].erase(i);
	i=find(connections[e.n2].begin(),connections[e.n2].end(),e.n1);
	connections[e.n2].erase(i);
}

void remove_lowest_edge(map<ID_T,list<ID_T> > &connections,map<edge,edge_data> &edges,set<ID_T> &x)
{	edge	min_edge;
	int val=numeric_limits<int>::max();
	for (set<ID_T>::iterator i=x.begin(); i!=x.end(); ++i)
	{	for (list<ID_T>::iterator j=connections[*i].begin(); j!=connections[*i].end(); ++j)
		{	edge	e(*i,*j);
			if (edges[e].count<val)
			{	min_edge=e;
				val=edges[e].count;
			}	
		}
	}
	remove_from_connections(connections,min_edge);
}

int main(int argc,char **argv)
{	string	fastafile;
	string	samfile;
	string	outfile;
	int	verbose=0;

	int k=0;

	while (-1!=(k=getopt(argc,argv,"hs:o:vf:")))
	{	switch (k)
		{	case 'h' :
				usage(argv[0]);
				return 0;
			case 's' :
				samfile=string(optarg);
				break;
			case 'f' :
				fastafile=string(optarg);
				break;
			case 'o' :
				outfile=string(optarg);
				break;
			case 'v' :
				verbose++;
				break;
			default :
				usage(argv[0]);
				return 7;
		}
	}

	if (samfile=="")
	{	cerr << "no graph input" << endl;
		return 1;
	}	

	if (fastafile=="")
	{	cerr << "no fasta input" << endl;
		return 1;
	}

	ifstream sam(samfile.c_str());
	if (!sam || !sam.is_open())
	{	cerr << "can't open '" << samfile << "'" << endl;
		return 2;
	}	

#if ID_TYPE == 1
	map<string,size_t> contig_sizes;
	get_contig_sizes(fastafile,contig_sizes);
#else
	map<string,size_t> h;
	get_contig_sizes(fastafile,h);
	map<int,size_t> contig_sizes;
	for (map<string,size_t>::iterator i=h.begin(); i!=h.end(); ++i)
	{	contig_sizes[lexical_cast<int>(i->first)]=i->second;
	}
	h.clear();
#endif	

	cout << contig_sizes.size() << " contigs in '" << fastafile << "'" << endl;

	map<edge,edge_data>	edges;
	ID_T	n1,n2;
	int count=0;
	int dist=0;
	while (true)
	{	string 	line;
		getline(sam,line);
		if (sam.eof())
		{	break;
		}
		istringstream is(line);
		string	s;
		is >> s;
		if (s=="edge:")
		{	is >> n1 >> n2 >> count;
			if (count==0)
			{	cerr << "count is 0" << endl;
				return 3;
			}	
		}
		else if (s=="dist")
		{	is >> dist;
		//	if (dist/count>=2)
		//	{
					edge e(n1,n2);
				if (edges.find(e)!=edges.end())
				{	cerr << "edge " << n1 << '-' << n2 << " already exists" << endl;
					return 4;
				}	
				//edges[e]=edge_data(count);
				edges.insert(pair<edge,edge_data>(e,edge_data(count)));
		//	}
		}	
		else if (s=="aligned")
		{	// silently ignored
		}
		else if (s=="1direct" || s=="2direct")
		{	// silently ignored
			int dodo;
			is >> dodo;
		}
		else
		{	cerr << "unknown tag '" << s << "'" << endl;
			return 5;
		}	
	}	
	cout << "edges read : " << edges.size() << endl;

	// check that we have a contig size for all used contigs
	bool no_contig_size=false;
	for (map<edge,edge_data>::iterator i=edges.begin(); i!=edges.end(); ++i)
	{	if (contig_sizes.find(i->first.n1)==contig_sizes.end())
		{	cerr << "no size for contig '" << i->first.n1 << "'" << endl;
			no_contig_size=true;
		}
		if (contig_sizes.find(i->first.n2)==contig_sizes.end())
		{	cerr << "no size for contig '" << i->first.n2 << "'" << endl;
			no_contig_size=true;
		}
	}
	if (no_contig_size)
	{	return 5;
	}

	// converting data structure
	map<ID_T,list<ID_T> >	connections;
	for (map<edge,edge_data>::iterator i=edges.begin(); i!=edges.end(); ++i)
	{	connections[i->first.n1].push_back(i->first.n2);
		connections[i->first.n2].push_back(i->first.n1);
	}

	// reducing until connections <=2 by :
	//	if contig size > any contig size of connected -> remove all connections
	//	if 3 or more connected contigs larger than self -> remove all connections
	//	if 2 connected contigs larger -> remove all smaller connections
	// 	if 1 connection -> keep it
	// 	generally : keep 3 largest including self

	unsigned int max_connections=0;
	for (map<ID_T,list<ID_T> >::iterator i=connections.begin(); i!=connections.end(); ++i)
	{	if (i->second.size()>max_connections)
		{       max_connections=i->second.size();
		}
	}	
	if (max_connections>2)
	{	do
		{	set<edge>	all_to_remove;
			for (map<ID_T,list<ID_T> >::iterator i=connections.begin(); i!=connections.end(); ++i)
			{	if (i->second.size()==max_connections)
				{	unsigned int self_size=contig_sizes[i->first];	
					list<ID_T>	larger;	// really >=
					list<ID_T>	smaller;
					for (list<ID_T>::iterator j=i->second.begin(); j!=i->second.end(); ++j)
					{	if (contig_sizes[*j]>=self_size)
						{	larger.push_back(*j);
						}
						else
						{	smaller.push_back(*j);
						}
					}
					if (smaller.size()==0 || larger.size()>=3)
					{	// remove all connections
						for (list<ID_T>::iterator j=i->second.begin(); j!=i->second.end(); ++j)
						{	all_to_remove.insert(edge(i->first,*j));

						}
					}
					if (larger.size()==2)
					{	// remove all smaller
						for (list<ID_T>::iterator j=smaller.begin(); j!=smaller.end(); ++j)
						{	all_to_remove.insert(edge(i->first,*j));
						}
					}	
				}	
			}
			for (set<edge>::iterator i=all_to_remove.begin(); i!=all_to_remove.end(); ++i)
			{	remove_from_connections(connections,*i);
			}
			max_connections--;
			cout << "max connections " << max_connections << endl;
		} while (max_connections>2);
	}	
	cout << "edges after removing connections > 2 : " << get_count_edges(connections) << endl;

	// circular detection
	int count_loops=0;
	for (map<ID_T,list<ID_T> >::iterator i=connections.begin(); i!=connections.end(); ++i)
	{	if (i->second.size()<2)
		{	continue;
		}
		set<ID_T>	visited;
		set<edge>	done;
		stack<edge>	todo;

		visited.insert(i->first);
		todo.push(edge(i->first,i->second.front()));
		todo.push(edge(i->first,i->second.back()));

		while (todo.size())
		{
			edge	e(todo.top());
			todo.pop();
			done.insert(e);
			ID_T from,to;
			if (visited.find(e.n1)!=visited.end())
			{	from=e.n1;
				to=e.n2;
			}
			else
			{	from=e.n2;
				to=e.n1;
			}
			if (visited.find(to)!=visited.end())
			{	cout << "loop detected : " << e << endl;
				count_loops++;
				remove_lowest_edge(connections,edges,visited);
				break;
			}
			else
			{	visited.insert(to);
				edge ne(to,connections[to].front());
				if (done.find(ne)==done.end())
				{	todo.push(ne);
				}
				if (connections[to].size()==2)
				{	ne=edge(to,connections[to].back());
					if (done.find(ne)==done.end())
					{       todo.push(ne);
					}
				}	
			}
		}
	}
	cout << count_loops << " loops broken" << endl;
	cout << "edges after breaking loops : " << get_count_edges(connections) << endl;

	// output
	set<ID_T>	visited;
	if (outfile!="")
	{	cout << "writing '" << outfile << "'" << endl;
		ofstream out(outfile.c_str());
		if (!out || !out.is_open())
		{	cerr << "can't open '" << outfile << "' for writing" << endl;
			return 5;
		}	
		
		for (map<ID_T,list<ID_T> >::iterator i=connections.begin(); i!=connections.end(); ++i)
		{	if (i->second.size()==0)
			{	out << i->first << endl;
			}
			else if (i->second.size()==1)
			{	if (visited.find(i->first)!=visited.end())
				{	continue;
				}

				out << i->first;
				visited.insert(i->first);
				ID_T c;
				c=i->second.front();
				while (c!=ID_T())
				{	out << '\t' << c;
					visited.insert(c);
					ID_T	n;
					n=connections[c].front();
					if (visited.find(n)==visited.end())
					{	c=n;
						continue;
					}	
					if (connections[c].size()==2)
					{	n=connections[c].back();
						if (visited.find(n)==visited.end())
						{	c=n;
							continue;
						}
					}
					c=ID_T();
				}	
				out << endl;
			}
		}
	}

	return 0;
}


void get_contig_sizes(const string &infile,map<string,size_t> &result)
{   int     in=open(infile.c_str(),O_RDONLY);
    size_t  input_size=lseek(in,0,SEEK_END);
    lseek(in,0,SEEK_SET);
    const size_t    part_size=128*1024*1024;
    char    *buffer=new char[part_size];
    size_t  pos=0;
    string  c_name;
    size_t  c_size=0;
    bool    in_header=false;
    while (pos<input_size)
    {       ssize_t readed=read(in,buffer,part_size);
        for (ssize_t i=0; i<readed; ++i)
        {   if (in_header)
            {       if (buffer[i]=='\n')
                {       in_header=false;
                    c_size=0;
                }
                else
                {   c_name+=buffer[i];
                }
            }
            else
            {   if (buffer[i]=='>')
                {       if (c_name!="")
                    {   result[c_name]=c_size;
                        c_name="";
                    }
                    in_header++;
                }
                else
                {   if (buffer[i]!='\n' && buffer[i]!='\t' && buffer[i]!=' ')
                    {   c_size++;
                    }
                }
            }
        }
        pos+=readed;
    }
    if (c_name!="")
    {   result[c_name]=c_size;
    }
}

