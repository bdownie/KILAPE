#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>

#include <getopt.h>

#define RIGHT 1
#define LEFT 0


using namespace std;

bool debug = 0;

class edge {
		unsigned int vertex1;
		unsigned int vertex2;
		unsigned int count;
		int distance;
		unsigned int direct1;
		unsigned int direct2;
	public:
		void v1(unsigned int x) { vertex1 = x; }
		void v2(unsigned int x) { vertex2 = x; }	
		void c(unsigned int x) { count = x; }
		void d1(unsigned int x) { direct1 = x; }	
		void d2(unsigned int x) { direct2 = x; }	
		void dist(int x) { distance = x; }

		unsigned int v1(void) { return vertex1; }
		unsigned int v2(void) { return vertex2; }
		unsigned int c(void) { return count; }
		unsigned int d1(void) { return direct1; }
		unsigned int d2(void) { return direct2; }
		int dist(void) { return distance; }

		bool is_edge(edge);
		bool is_edge(unsigned int, unsigned int);

		void print_formatted(void);
		void print_orphan(void);
		void init(void);
		void assign(unsigned int, unsigned int, unsigned int, int, unsigned int, unsigned int); 
		void debug(void);
};



void edge::debug(void) {
	cout << vertex1 << "\t";
	cout << vertex2 << "\t";
	cout << count << "\t";
	cout << distance << "\t";
	cout << direct1 << "\t";
	cout << direct2 << endl;
}

void edge::assign(unsigned int a, unsigned int b , unsigned int c, int d, unsigned int e, unsigned int f) {
	vertex1 = a;
	vertex2 = b;
	count = c;
	distance = d;
	direct1 = e;
	direct2 = f;
}

void edge::init(void) { 
	vertex1 = 0;
	vertex2 = 0;
	count = 0;
	distance = 0;
	direct1 = -1;
	direct2 = -1;
}

void edge::print_orphan(void) { 
//	string return_string;
//	stringstream tmp;
	ofstream orphan;
	orphan.open("orphan.graph", ios::app);

	orphan << "edge: " << vertex1 << " " << vertex2 << " " << count << "\n";
	orphan << "dist " << distance << "\n";
	orphan << "1direct " << direct1 << "\n";
	orphan << "2direct " << direct2 << "\n";

	orphan.close();
//	return return_string;
}

void edge::print_formatted(void) { 
	cout << "edge: " << vertex1 << " " << vertex2 << " " << count << endl;
	cout << "dist " << distance << endl;
	cout << "1direct " << direct1 << endl;
	cout << "2direct " << direct2 << endl;
}

class vertex {
		unsigned int name_private;
		int size_private;

		// index is the destination vertex
		map <unsigned int, edge> edges;
		map <unsigned int,edge>::iterator e_it;
		
	public:
		unsigned int name(void) { return name_private; }
		void name(unsigned int x) { name_private = x; }

		void size(int x) { size_private = x; }
		int size(void) { return size_private; }

		bool edge_direction(edge);
		bool is_vertex(unsigned int);

		unsigned int get_dest(edge);
		unsigned int get_dest(unsigned int, unsigned int);


		edge get_edge(unsigned int x) { return edges[x]; }
		void insert_edge(edge);
		void delete_edge(edge);
		void delete_edge(unsigned int, unsigned int);

		vector <edge> get_all_edges(void);
		map <int, edge> get_directional_edges(bool);
		void debug(void);
};


map<int, edge> vertex::get_directional_edges(bool b) { 
	// b = 1 for right
	// b = 0 for left

	// return index is distance between stuff.

	//vector <edge> return_vector;

	map <int, edge> return_map;
	map <int, edge> rm_it;

	for (e_it = edges.begin(); e_it != edges.end(); e_it++) { 
		bool bool_d1 = ((*e_it).second.d1() != 0);
		bool bool_d2 = ((*e_it).second.d2() != 0);
		int distance = (*e_it).second.dist();
		if ((*e_it).second.v1() == name_private) {
			if (bool_d1 == b) {
				return_map.insert(pair<int,edge>(distance,(*e_it).second));
			}
		}
		else if (bool_d2 == b) {
			return_map.insert(pair<int,edge>(distance,(*e_it).second));
		}
	}

	return return_map;
}

void vertex::debug(void) { 
	for (e_it = edges.begin(); e_it != edges.end(); e_it++) { 
		(*e_it).second.debug();
	}
}


unsigned int vertex::get_dest(unsigned int x, unsigned int y) { 
	unsigned int dest = 0;
	if (name_private == x) { dest = y; }
	else if (name_private == y) { dest = x; }
	e_it = edges.find(dest);
	if (e_it == edges.end()) { dest = 0; }
	//else { cout << "Edge doesn't belong to vertex. " << x << " "<< y << endl; }

	return dest;
}

unsigned int vertex::get_dest(edge e) { 
	unsigned int dest = 0;
	if (name_private == e.v1()) { 
		dest = e.v2();
	}
	else if (name_private == e.v2()) { 
		dest = e.v1();
	}
	//else { cout << "Edge doesn't belong to vertex: " <<name_private << " " << e.v1() << " " << e.v2()<< endl; }

	return dest;
}

void vertex::delete_edge(edge e) {
	unsigned int dest = get_dest(e);
	edges.erase(dest);
}

void vertex::delete_edge(unsigned int x, unsigned int y) {
	unsigned int dest = get_dest(x,y);
	edges.erase(dest);
}

void vertex::insert_edge(edge e) { 
	unsigned int dest = get_dest(e);
	int distance = e.dist();

	if (edges.size() > 0) { 
		for (e_it = edges.begin(); e_it != edges.end(); e_it++) { 
			if ((*e_it).second.is_edge(e)) { return; }
			if ((*e_it).second.dist() == distance) { 
				e.dist(++distance);
			}
		}
	}
	edges.insert(pair<unsigned int, edge>(dest,e));
//	if (e_it == edges.end()) { cerr << "Failed to insert edge" << endl; }
}

vector <edge> vertex::get_all_edges(void) {
	vector <edge> return_vector;
	for (e_it = edges.begin(); e_it != edges.end(); e_it++) { 
		return_vector.push_back((*e_it).second);
	}
	return return_vector;
}

unsigned int orphan_count = 0;


void printUsage();
map<unsigned int, vertex> readGraph(char *, char *, unsigned int);
vertex get_connecting_vertex(map <unsigned int, vertex> *v, unsigned int source, edge e);
bool create_new_edge(map <unsigned int, vertex> *, unsigned int , unsigned int, unsigned int , unsigned int , int);
void add_edge_to_vertex( map<unsigned int, vertex> *, unsigned int ,   edge);
int delete_edge(vector <edge> *, edge , string );
void print_new_graph_file(map <unsigned int, vertex> *v);
void simplify_graph(map <unsigned int, vertex> *v,unsigned int);
edge delete_outer_edge(map <unsigned int, vertex> *, map <int,edge> *, unsigned int, unsigned int );
void create_inner_edge(map <unsigned int, vertex> *, map <int, edge> *, edge );
unsigned int get_threshold_size(map<int,edge> , unsigned int );
void print_scaff_file(map<unsigned int, vertex> *v, char *);
unsigned int get_dest (unsigned int source, edge e);
void trim_large_inserts(map <unsigned int, vertex> *v);
map <int, edge>  trim_large_inserts_edges(map <int, edge>, map <unsigned int, vertex> *v);

int main(int argc, char* argv[]) {
	char c;
	char *GraphFile = NULL;
	char *ContigFile =  NULL;
	char *ScaffFile =  NULL;
	//debug = 0;
	unsigned int min_pair = 2;

	while ((c = getopt (argc, argv, "dhg:f:p:s:")) != -1) {
		switch(c) {
			case 's':
				ScaffFile = optarg;
				break;
			case 'p':
				min_pair = atoi(optarg);
				break;
			case 'g':
				GraphFile = optarg;
				break;
			case 'f':
				ContigFile = optarg;
				break;
			case 'h':
				printUsage();
				return 0;
			case 'd':
				debug = 1;
				break; 
			default:
				printUsage();
				return 0;
		}
	}

	map <unsigned int, vertex> vertices;
	
	if ((GraphFile == NULL) || (ContigFile == NULL) ) { printUsage(); return 0; }
	vector <edge> all_edges;
	cerr << "Reading graph file" << endl;
	vertices = readGraph(GraphFile, ContigFile, min_pair);

/*
cout << "debug" << endl;
	vertices[262729].debug();
cout << "debug2" << endl;
	vertices[185241].debug();
	cout << endl;
	*/
	cerr << "Finished reading initial graph discrepencies" << endl;

	
	//trim_large_inserts(&vertices);
	//return 0;
	simplify_graph(&vertices, min_pair);
	print_new_graph_file(&vertices);
	cerr << endl << "Orphan Count: " << orphan_count << endl;

	if (ScaffFile != NULL) { 
		print_scaff_file(&vertices, ScaffFile);
	}
	

	cerr <<endl;

	return 0;
}

void print_scaff_file(map<unsigned int, vertex> *v, char *file) {
	map <unsigned int, vertex>::iterator v_it;
	vector<unsigned int> visited;

	v_it = (*v).end();

	visited.resize((*v_it).first, 0);

	ofstream scaffold;
	scaffold.open(file);

/*
cout << "debug" << endl;
	(*v)[262729].debug();
cout << "debug2" << endl;

	(*v)[185241].debug();
cout << endl;	
*/
	for (v_it = (*v).begin(); v_it != (*v).end(); v_it++) {
	//	scaffold << "Working on..." << (*v_it).first << endl;
		if (visited[(*v_it).first]) { continue; }

		vector <edge> v_edges = (*v_it).second.get_all_edges();
		
		if (v_edges.size() == 0) { 
			scaffold << (*v_it).first << endl; 
			visited[(*v_it).first] = 1;
		}
		else if (v_edges.size() == 1) { 


			unsigned int current_vertex = (*v_it).first;
			unsigned int next_vertex = get_dest(current_vertex, v_edges[0]);
			//(*v_it).second.debug();
			//v_edges[0].debug();
	
	
	//scaffold << current_vertex << "-" << next_vertex << " ";
			scaffold << (*v_it).first << "\t" << next_vertex;
		
			v_edges = (*v)[next_vertex].get_all_edges();
			visited[current_vertex] = 1;
			while ((visited[next_vertex] == 0) && (v_edges.size() > 1)) {
				current_vertex = next_vertex;
				next_vertex =  get_dest(current_vertex, v_edges[0]);
				if (visited[next_vertex]) { 
					next_vertex =  get_dest(current_vertex, v_edges[1]);
				}
				//cout << " " << current_vertex << "-" << next_vertex;
				//(*v)[current_vertex].debug();
	
				scaffold << "\t" << next_vertex;
				visited[current_vertex] = 1;
				v_edges = (*v)[next_vertex].get_all_edges();
			}
			scaffold << "\n";
		}
	}
	/*
	for (unsigned int i = 0; i < visited.size(); i++) { 
		if (visited[i]) { continue; }

		vector <edge> v_edges = (*v)[i].get_all_edges();
		unsigned int current_vertex = i;
		unsigned int next_vertex = get_dest(current_vertex, v_edges[0]);

	
		scaffold << i << "\t" << next_vertex;

		v_edges = (*v)[next_vertex].get_all_edges();
		visited[current_vertex] = 1;
		while ((visited[next_vertex] == 0) && (v_edges.size() > 1)) {
			current_vertex = next_vertex;
			next_vertex =  get_dest(current_vertex, v_edges[0]);
			if (visited[next_vertex]) { 
				next_vertex =  get_dest(current_vertex, v_edges[1]);
			}

			scaffold << "\t" << next_vertex;
			visited[current_vertex] = 1;
			v_edges = (*v)[next_vertex].get_all_edges();
		}
		scaffold << "\n";
	}
	*/
}

unsigned int get_dest (unsigned int source, edge e) {
	if (source == e.v1()) { 
		return e.v2();
	}
	else if (source == e.v2()) { 
		return e.v1();
	}
	else { 
		cout << "Bad edge attached to vertex " << source << endl;
		//e.debug();
		throw "error";
	}

	return 0;
}

void trim_large_inserts(map <unsigned int, vertex> *v) {
	map <unsigned int, vertex>::iterator v_it;


	for (v_it = (*v).begin(); v_it != (*v).end(); v_it++) {
		vector<edge> edges = (*v_it).second.get_all_edges();
		for (vector<edge>::iterator e_it = edges.begin(); e_it != edges.end(); e_it++) {
			unsigned int v1 = (*e_it).v1();
			unsigned int v2 = (*e_it).v2();
			int total_size = (*v)[v1].size() + (*v)[v2].size(); 
			total_size -= 200;
			if (total_size < (*e_it).dist()) {

				e_it->print_orphan();
				orphan_count++;
				//cout << "Deleting" << v1 << " " << v2 << " " << total_size << " " << (*e_it).dist() << endl;
				(*v_it).second.delete_edge(*e_it);
			}
			//if ((0 - total_size) > (*e_it).dist()) {
				//cout << "Deleting" << v1 << " " << v2 << " " << total_size << " " << (*e_it).dist() << endl;
			//	(*v_it).second.delete_edge(*e_it);
			//}
		}
	}
}

map <int, edge> trim_large_insert_edges(map <int, edge> edge_map, map <unsigned int, vertex> *v) { 
	for (map<int,edge>::iterator edge_map_it = edge_map.begin(); edge_map_it != edge_map.end(); edge_map_it++) { 
		unsigned int v1 = edge_map_it->second.v1();
		unsigned int v2 = edge_map_it->second.v2();
		int distance = edge_map_it->second.dist();
		int total_size = (*v)[v1].size() + (*v)[v2].size(); 

		if (total_size < distance) { 
			edge_map.erase(edge_map_it);
			(*v)[v1].delete_edge(edge_map_it->second);
			(*v)[v2].delete_edge(edge_map_it->second);
			edge_map_it->second.print_orphan();
			orphan_count++;
		}

	}
	return edge_map;
}



void simplify_graph(map <unsigned int, vertex> *v, unsigned int threshold) {
	map <unsigned int, vertex>::iterator v_it;
	map <unsigned int, vertex>::iterator rv_it;


	map <unsigned int, vector<map<unsigned int, vertex>::iterator> > tocheck;
/*
	for (v_it = (*v).begin(); v_it != (*v).end(); v_it++) { 
		map <unsigned int, vector<map<unsigned int, vertex>::iterator> >::iterator x;
		vector<edge> tmp_edges = (*v_it).get_all_edges();
		unsigned int count = tmp_edges.size();
		x = tocheck.find(count);
		if (x == tocheck.end()) {
			tocheck.insert(pair<unsigned int, vector<map<unsigned int, vertex>::iterator>
			*/

/*
	vector<vector<unsigned int> > count_sorted_vertices;
	vector<vector<unsigned int> >::iterator csv_it;
	unsigned int full_size = (*v).size();
	full_size++;
	//count_sorted_vertices.resize(full_size);

	for (v_it = (*v).begin(); v_it != (*v).end(); v_it++) { 
		vector<edge> tmp_edges = (*v_it).second.get_all_edges();
		if (tmp_edges.size() < 2) { continue; }
	 	unsigned int max_count = 0;
		for (unsigned int i = 0; i < tmp_edges.size(); i++) {
			if (tmp_edges.size() > max_count) { 
				max_count = tmp_edges.size();
			}
		}
		unsigned int name = (*v_it).second.name();
		if (count_sorted_vertices.size() <= max_count) {
			count_sorted_vertices.resize(max_count+1);
		}
		count_sorted_vertices[max_count].push_back(name);
	}
	*/

	cerr << "Ordering contigs by size" << endl;
	unsigned int full_size = (*v).size();
	full_size++;
	vector <vector<unsigned long> > dist_sorted_vertices;
	vector <vector<unsigned long> >::iterator dsv_it;
	dist_sorted_vertices.resize(full_size);
	for (v_it = v->begin(); v_it != v->end(); v_it++) { 
		unsigned long length = 0;

		length = v_it->second.size();
		if (length > dist_sorted_vertices.size()) { 
			dist_sorted_vertices.resize(length + 1);
		}

//cout << v_it->second.name() << " has length " << length << endl;
		dist_sorted_vertices[length].push_back(v_it->second.name());
	}
	cerr << "Done ordering contigs" << endl;



	unsigned int numVerts = 0;
	unsigned int lastSize = 1;
	cerr << endl << "Processing contigs: 0";
	//for (vector<unsigned int >::reverse_iterator i = dist_sorted_vertices.rbegin(); i != dist_sorted_vertices.rend(); i++) {
	for (long i = dist_sorted_vertices.size() - 1; i >= 0; i--) { 
	//for (long i = 0; i < dist_sorted_vertices.size(); i++) { 
		//if ((*i) == 0) { continue; }
		if (dist_sorted_vertices[i].empty()) { continue; }
	//	cout << endl << "Processing length of " << i << ". Set has size: " << dist_sorted_vertices[i].size();
		for (vector<unsigned long>::iterator j = dist_sorted_vertices[i].begin(); j != dist_sorted_vertices[i].end(); j++) { 
	//	cout << endl << "Processing contig: " << *j << " with length " << i;
		//v_it = (*v).find(*i);
		v_it = v->find(*j);

		numVerts++;

		if ((numVerts % 1000) == 0) { 
			for (unsigned long k = 0; k < lastSize; k++) {
				cerr << "\b";
			}
			stringstream Output;
			string strOutput;
			Output << numVerts;
			Output >> strOutput;
			lastSize = strOutput.size();
			cerr << numVerts;
		}
		unsigned int vertex_name = (*v_it).first;

		bool trim_all_edges = 0;
		map <int,edge> edges_to_trim = (*v_it).second.get_directional_edges(LEFT);
		if (edges_to_trim.size() > 1){
	//		edges_to_trim = trim_large_inserts_edges(edges_to_trim, v);
	//		run_again = 1; 
	/*
			if (edges_to_trim.size() > 2) { 
				for (map<int,edge>::iterator i = edges_to_trim.begin(); i != edges_to_trim.end(); i++) { 
					//(*v_it).second.delete_edge((*i).second);
					unsigned int a, b;
					a = (*i).second.v1();
					b = (*i).second.v2();
					(*v)[a].delete_edge((*i).second);
					(*v)[b].delete_edge((*i).second);
					edges_to_trim.erase(i);
				}
			}
			*/
			while (edges_to_trim.size() > 1){
/*

						*/
				edge new_edge = delete_outer_edge(v,&edges_to_trim, vertex_name, threshold);
				if ((new_edge.v1() != 0) && (new_edge.v2() != 0) && (!trim_all_edges)) { 
					create_inner_edge(v, &edges_to_trim,new_edge);
				}
			}
		}
			trim_all_edges = 0;
			/*
			while (edges_to_trim.size() > 1){
				edge new_edge = delete_outer_edge(v,&edges_to_trim, vertex_name, threshold);
//					if (threshold == 0) {
//						if (new_edge.dist() < 20) { new_edge.dist(20); }
//					}
				create_inner_edge(v, &edges_to_trim,new_edge);
			}
			*/
	
/*
			if (edges_to_trim.size() > 1) { 
				unsigned int threshold = get_threshold_size(edges_to_trim, 5);
				for (map<int,edge>::iterator i = edges_to_trim.begin(); i != edges_to_trim.end(); i++) { 
					if ((*i).second.c() <= threshold) {
						if (((*i).second.dist() == 2) && ((*i).second.c() == 1)) {
						//	(*v_it).second.delete_edge((*i).second);
						unsigned int a, b;
						a = (*i).second.v1();
						b = (*i).second.v2();
						(*v)[a].delete_edge((*i).second);
						(*v)[b].delete_edge((*i).second);
						}
					}
				}
				edges_to_trim = (*v_it).second.get_directional_edges(RIGHT);
				if (edges_to_trim.size() > 1) { 
					for (map<int,edge>::iterator i = edges_to_trim.begin(); i != edges_to_trim.end(); i++) { 
						//(*v_it).second.delete_edge((*i).second);
						unsigned int a, b;
						a = (*i).second.v1();
						b = (*i).second.v2();
						(*v)[a].delete_edge((*i).second);
						(*v)[b].delete_edge((*i).second);
	
					}
				}
				edges_to_trim = (*v_it).second.get_directional_edges(RIGHT);
			}

*/
			edges_to_trim = (*v_it).second.get_directional_edges(RIGHT);
			if (edges_to_trim.size() > 1){
	//			edges_to_trim = trim_large_inserts_edges(edges_to_trim, v);
	//		run_again = 1; 
	/*
				if (edges_to_trim.size() > 2) { 
					for (map<int,edge>::iterator i = edges_to_trim.begin(); i != edges_to_trim.end(); i++) { 
						//(*v_it).second.delete_edge((*i).second);
						unsigned int a, b;
						a = (*i).second.v1();
						b = (*i).second.v2();
						(*v)[a].delete_edge((*i).second);
						(*v)[b].delete_edge((*i).second);
						edges_to_trim.erase(i);
					}
				}
				*/
				while (edges_to_trim.size() > 1){
				/*
				while (edges_to_trim.size() > 0){
						for (map<int,edge>::iterator i = edges_to_trim.begin(); i != edges_to_trim.end(); i++) { 
							//(*v_it).second.delete_edge((*i).second);
							unsigned int a, b;
							a = (*i).second.v1();
							b = (*i).second.v2();
							(*v)[a].delete_edge((*i).second);
							(*v)[b].delete_edge((*i).second);
							edges_to_trim.erase(i);
		
						}
						*/
					edge new_edge = delete_outer_edge(v,&edges_to_trim, vertex_name, threshold);
					if ((new_edge.v1() != 0) && (new_edge.v2() != 0) && (!trim_all_edges)) { 
						create_inner_edge(v, &edges_to_trim,new_edge);
					}
				}
			}
		}
		}
	//	if (run_again) { 
	//		cout << endl << "Running again..." << endl; 
	//	}
}

unsigned int get_threshold_size(map<int,edge> ie, unsigned int divisor) {
	unsigned int result = 0;
	for (map<int,edge>::iterator i = ie.begin(); i != ie.end(); i++) { 
		unsigned int count_tmp = (*i).second.c();
		if (count_tmp > result) { 
			result = count_tmp;
		}
	}

	result = (unsigned int)(result/divisor);
	result++;

	return result;
}

void create_inner_edge(map <unsigned int, vertex> *v, map <int, edge> *ett, edge e) { 
	map <int,edge>::iterator ett_it;
	map <int,vertex>::iterator v_it;
	/*
	if (e.v1() == e.v2()) {
		return;
	}
	*/

	for (ett_it = (*ett).begin(); ett_it != (*ett).end(); ett_it++) { 
		if ((*ett_it).second.is_edge(e)) {
			return;
		}
	}

	ett_it = (*ett).find(e.dist());
	while (ett_it != (*ett).end()) {
		int dist = e.dist() + 1;
		e.dist(dist);
		ett_it = (*ett).find(e.dist());
	}

	
//	cout << "Creating at vector" << endl;
//e.debug()
	//(*ett).insert(pair<int, edge>(e.dist(),e));
	(*v)[e.v1()].insert_edge(e);
	(*v)[e.v2()].insert_edge(e);
}

edge delete_outer_edge(map <unsigned int, vertex> *v, map <int,edge> *ett, unsigned int source_name, unsigned int threshold) {
	map <int, edge>::iterator outside_edge;
	map <int, edge>::iterator inside_edge;


	edge trash_edge;
	trash_edge.assign(0,0,0,0,0,0);
	edge return_edge = trash_edge;

	outside_edge = (*ett).end();
	outside_edge--;
	inside_edge = outside_edge;
	inside_edge--;

	unsigned int outside_name, inside_name;
	int outside_edge_dir, inside_edge_dir;

	


//inside_edge->second.debug();
//outside_edge->second.debug();
	if (source_name == (*outside_edge).second.v1()) {
		outside_name = (*outside_edge).second.v2();
		outside_edge_dir = (*v)[outside_name].edge_direction((*outside_edge).second);
//	cout << "here1 " << source_name << " " << outside_name << " " << outside_edge_dir << endl;
//		outside_edge->second.debug();
	}
	else { 
		outside_name = (*outside_edge).second.v1();
		outside_edge_dir = (*v)[outside_name].edge_direction((*outside_edge).second);
//	cout << "here2 " << source_name << " " << outside_name << " " << outside_edge_dir << endl;
//		outside_edge->second.debug();
	}
//	cout << source_name << " " << outside_name << " " << outside_edge_dir << endl;
//	(*outside_edge).second.debug();

	if (source_name == (*inside_edge).second.v1()) {
		inside_name = (*inside_edge).second.v2();
		inside_edge_dir = (*v)[inside_name].edge_direction((*inside_edge).second);
	}
	else { 
		inside_name = (*inside_edge).second.v1();
		inside_edge_dir = (*v)[inside_name].edge_direction((*inside_edge).second);
	}
	//cout << outside_name << " " << outside_edge_dir << endl;
	map <unsigned int, vertex>::iterator v_it;
	v_it = (*v).find(inside_name);
	if (v_it != (*v).end()) {
		if ((*v_it).second.get_dest(inside_name, outside_name)) {
	//		cout << "Current Vertex : " << source_name << endl;
	//		cout << "Edge already exists: " << inside_name << " " << outside_name << endl;
			(*v)[outside_name].delete_edge((*outside_edge).second);
			(*v)[source_name].delete_edge((*outside_edge).second);
			(*ett).erase(outside_edge);
			//return (*inside_edge).second;
			return trash_edge;
		}
	}

//cout <<source_name << " " << inside_name << " " << inside_edge_dir << endl;
//cout << outside_name << " " << outside_edge_dir << endl;

	if ((((*v)[outside_name].size()) + ((*v)[source_name].size())) < outside_edge->second.dist()) {
		outside_edge->second.print_orphan();
		orphan_count++;
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return trash_edge;
	}
	else if ((((*v)[inside_name].size()) + ((*v)[source_name].size())) < inside_edge->second.dist()) {
		inside_edge->second.print_orphan();
		orphan_count++;
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		return (*outside_edge).second;
	}

	// Eliminate edges which have connections outside of threshold
	if ((*outside_edge).second.c() < threshold) {
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return (*inside_edge).second;
	}
	if ((*inside_edge).second.c() < threshold) {
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		return (*outside_edge).second;
	}

	if ((*v)[inside_name].size() > (outside_edge->second.dist() * 2)) { 
		if (inside_edge->second.c() >= outside_edge->second.c()) { 
			(*v)[outside_name].delete_edge((*outside_edge).second);
			(*v)[source_name].delete_edge((*outside_edge).second);
			(*ett).erase(outside_edge);
			return (*inside_edge).second;
		}
		else { 
			(*v)[inside_name].delete_edge((*inside_edge).second);
			(*v)[source_name].delete_edge((*inside_edge).second);
			(*ett).erase(inside_edge);
			return (*outside_edge).second;
		}
	}


	bool outside_edge_dir_bool = !(!(outside_edge_dir));
	bool inside_edge_dir_bool = !(inside_edge_dir);

	map <int, edge> outside_edges = (*v)[outside_name].get_directional_edges(outside_edge_dir_bool);
	map<int, edge> inside_edges = (*v)[inside_name].get_directional_edges(inside_edge_dir_bool);


	int new_edge_dist = (*outside_edge).second.dist() - ((*v)[inside_name].size() + (*inside_edge).second.dist());
	int new_edge_count = (int)(((*outside_edge).second.c() + (*inside_edge).second.c())/2);

	


	int min_contig_size = 0 - new_edge_dist;

	if ((min_contig_size > (*v)[inside_name].size())
						|| 
	    (min_contig_size > (*v)[outside_name].size())
						|| 
	    (min_contig_size > (*v)[source_name].size())) { 

		if ((outside_edge->second.c()) > (inside_edge->second.c() * 2)) {
			inside_edge->second.print_orphan();
			orphan_count++;
			(*v)[inside_name].delete_edge((*inside_edge).second);
			(*v)[source_name].delete_edge((*inside_edge).second);
			(*ett).erase(inside_edge);
			return (*outside_edge).second;
		}
		else {
			outside_edge->second.print_orphan();
			orphan_count++;
			(*v)[outside_name].delete_edge((*outside_edge).second);
			(*v)[source_name].delete_edge((*outside_edge).second);
			(*ett).erase(outside_edge);
			return trash_edge;
		}
	}

/*
	if ((inside_edge->second.c()) > (3 * outside_edge->second.c())) {
		outside_edge->second.print_orphan();
		orphan_count++;
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return trash_edge;
	}
	else if ((outside_edge->second.c()) > (inside_edge->second.c() * 3)) {
		inside_edge->second.print_orphan();
		orphan_count++;
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		return (*outside_edge).second;
	}
	*/

	vector<edge> edges_connected_to_inside;
	// Make sure distant contig is only connected to source contig
	// in the proper direction.
	map <int, edge> outside_edge_dir_edges;

//cout << endl << "Debug: " << endl;
//outside_edge->second.debug();
//inside_edge->second.debug();
//return_edge.debug();
	if (new_edge_count >= 1) { 
	//	cout << source_name << " " << inside_name << " " << outside_name << " " << inside_edge_dir << " " << outside_edge_dir << endl;
		int new_inside_edge_dir = 0;
		if (inside_edge_dir == 0) { new_inside_edge_dir = new_edge_count; }
		if (outside_edge_dir != 0) { outside_edge_dir = new_edge_count; }
	
		if (inside_name < outside_name) { 
			return_edge.assign(inside_name, outside_name, new_edge_count, new_edge_dist, new_inside_edge_dir, outside_edge_dir);
		}
		else {
			return_edge.assign(outside_name, inside_name, new_edge_count, new_edge_dist, outside_edge_dir, new_inside_edge_dir);
		}
		// Make sure the middle contig is only connected to source contig
		edges_connected_to_inside = (*v)[inside_name].get_all_edges();
		// Make sure distant contig is only connected to source contig
		// in the proper direction.
		outside_edge_dir_edges = (*v)[outside_name].get_directional_edges(outside_edge_dir_bool);
	}

	//cout << "Sizes: " << edges_connected_to_inside.size() << " " << outside_edge_dir_edges.size() << endl;
//	if ((edges_connected_to_inside.size() <= 1) && (outside_edge_dir_edges.size() <= 1)) {

	outside_edge->second.print_orphan();
	orphan_count++;
	(*v)[outside_name].delete_edge((*outside_edge).second);
	(*v)[source_name].delete_edge((*outside_edge).second);
	(*ett).erase(outside_edge);
		return return_edge;
//	}



/*	

	else { 
		inside_edge->second.print_orphan();
		orphan_count++;
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		outside_edge->second.print_orphan();
		orphan_count++;
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		
		return trash_edge;
	//	return (*inside_edge).second;
	}
	
	

	if (new_edge_dist < 0) {
		if ((*v)[outside_name].size() > (*v)[inside_name].size()) {
			inside_edge->second.print_orphan();
			orphan_count++;
			(*v)[inside_name].delete_edge((*inside_edge).second);
			(*v)[source_name].delete_edge((*inside_edge).second);
			(*ett).erase(inside_edge);
			return (*outside_edge).second;
		}
		else { 
			outside_edge->second.print_orphan();
			orphan_count++;
			(*v)[outside_name].delete_edge((*outside_edge).second);
			(*v)[source_name].delete_edge((*outside_edge).second);
			(*ett).erase(outside_edge);
			return (*inside_edge).second;
		}
	}



	if ((*inside_edge).second.c() >= ((*outside_edge).second.c() * 2)) {
		outside_edge->second.print_orphan();
		orphan_count++;
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return (*inside_edge).second;
	}
	if ((*outside_edge).second.c() >= ((*inside_edge).second.c() * 2)) {
		inside_edge->second.print_orphan();
		orphan_count++;
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		return (*outside_edge).second;
	}
	*/

	/*
	*/

	/*
	if ((outside_edges.size() > 1) || (inside_edges.size() > 0)) {
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return (*inside_edge).second;
	}
	*/


/*
	if ((*outside_edge).second.dist() == 2) { 
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return (*inside_edge).second;
	}
	else if ((*inside_edge).second.dist() == 2) { 
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		return (*outside_edge).second;
	}
	*/
	/*
	if ((*outside_edge).second.c() > (5 * (*inside_edge).second.c())) {
		(*v)[inside_name].delete_edge((*inside_edge).second);
		(*v)[source_name].delete_edge((*inside_edge).second);
		(*ett).erase(inside_edge);
		//cout << "deleting " << (*inside_edge).second.v1() << " " << (*inside_edge).second.v2() << endl;
		return (*outside_edge).second;
	}
	*/
/*
	inside_edge_dir = !inside_edge_dir;

	if ((new_edge_dist < (0 - (*outside_edge).second.dist())) || (new_edge_count == 1)) { 
		if ((*outside_edge).second.c() < (*inside_edge).second.c()) { 
			inside_edge->second.print_orphan();
			orphan_count++;
			(*v)[inside_name].delete_edge((*inside_edge).second);
			(*v)[source_name].delete_edge((*inside_edge).second);
			(*ett).erase(inside_edge);
			return (*outside_edge).second;
		}
		else { 
			outside_edge->second.print_orphan();
			orphan_count++;
			(*v)[outside_name].delete_edge((*outside_edge).second);
			(*v)[source_name].delete_edge((*outside_edge).second);
			(*ett).erase(outside_edge);
			return (*inside_edge).second;
		}
	}
	else if (new_edge_count < 1) { 
		outside_edge->second.print_orphan();
		orphan_count++;
		(*v)[outside_name].delete_edge((*outside_edge).second);
		(*v)[source_name].delete_edge((*outside_edge).second);
		(*ett).erase(outside_edge);
		return (*inside_edge).second;
	}

*/
	if (inside_edge_dir != 0) { inside_edge_dir = new_edge_count; }
	if (outside_edge_dir != 0) { outside_edge_dir = new_edge_count; }

	if (inside_name < outside_name) { 
		return_edge.assign(inside_name, outside_name, new_edge_count, new_edge_dist, inside_edge_dir, outside_edge_dir);
	}
	else {
		return_edge.assign(outside_name, inside_name, new_edge_count, new_edge_dist, outside_edge_dir, inside_edge_dir);
	}

	outside_edge->second.print_orphan();
	orphan_count++;
	(*v)[outside_name].delete_edge((*outside_edge).second);
	(*v)[source_name].delete_edge((*outside_edge).second);
	(*ett).erase(outside_edge);

	return return_edge;
}

void print_new_graph_file(map <unsigned int, vertex> *v) {

	map <unsigned int, vertex>::iterator v_it;
	map <unsigned int, vertex>::iterator rv_it;

//cout << "graph file size: " << (*v).size() << endl;
	for (v_it = (*v).begin(); v_it != (*v).end(); v_it++) { 
		vector <edge> tmp_edges;

		vector <edge>::iterator ve_i;
		vector <edge>::iterator te_i;

		tmp_edges = (*v_it).second.get_all_edges();	
//		cout << "Tmp: " << tmp_edges.size() << endl;
		for (te_i = tmp_edges.begin(); te_i != tmp_edges.end(); te_i++) { 
			unsigned int v1 = (*te_i).v1();
			unsigned int v2 = (*te_i).v2();
			if ((*v_it).second.name() == v1) {
				rv_it = (*v).find(v2);
			}
			else if ((*v_it).second.name() == v2) {
				rv_it = (*v).find(v1);
			}
			else {
				cerr << "Unknown edge found ";
				//(*v_it).second.debug();
				cerr << endl;
				//return -1;
			}
			(*te_i).print_formatted();
			(*rv_it).second.delete_edge(v1, v2);
			//(*v_it).second.delete_edge(v1, v2);
		}
	}
}


map<unsigned int, vertex> readGraph(char *Graph, char *Contig, unsigned int minimum_pair) {
	ifstream fhGraph, fhContig;
	fhGraph.open(Graph);
	fhContig.open(Contig);

	unsigned int ep = 0;
//	vector<unsigned int> contig_sizes;

	map <unsigned int, vertex> verts;
	map <unsigned int, vertex>::iterator v_it;
	map <unsigned int, vertex>::iterator v2_it;

	if (!fhGraph.is_open()) {
		cerr  << "Couldn't open file: " << Graph<< endl;
		return verts;
	}
	if (!fhContig.is_open()) {
		cerr  << "Couldn't open file: " << Contig<< endl;
		return verts;
	}

	edge e;
	e.init();


	unsigned int output_size = 0;
	unsigned int acc = 0;
	while (fhContig.good()) { 
		string Line;
		//unsigned int output_size;
		getline(fhContig, Line);
		if (!fhContig.good()) { break; }

		if (Line[0] == '>') { 
			if (output_size) { 
				vertex tmp_vertex;
				tmp_vertex.size(output_size);
				tmp_vertex.name(acc);
				verts.insert(pair<unsigned int, vertex>(acc,tmp_vertex));
			}

			Line.erase(0,1);
				//acc = Line;
			acc = atoi(Line.c_str());
			//contig_sizes.push_back(output_size);
			output_size = 0;
		}
		else {
			output_size += Line.size();
		}
	}
	if (output_size) { 
		vertex tmp_vertex;
		tmp_vertex.size(output_size);
		tmp_vertex.name(acc);
		verts.insert(pair<unsigned int, vertex>(acc,tmp_vertex));
	}


	cerr << "Finished reading contigs" << endl;

	ofstream orphan;
	orphan.open("orphan.graph",ios::out);
	string orphan_string;
	while (fhGraph.good()) { 
		stringstream tmp;
		string Line;
		string edgeline;

		string param;
		unsigned int v1, v2;
		unsigned int count;
		unsigned int dist;
		unsigned int d1, d2;



		getline(fhGraph, Line);
		if (!fhGraph.good()) { break; }

		tmp << Line;
		tmp >> param;
		orphan_string.append(Line);
		orphan_string.append("\n");

		if (param == "edge:") { 
			edgeline = Line;
			tmp >> v1;
			tmp >> v2;
			tmp >> count;
			if (v2 < v1) { 
				unsigned int tmp2 = v2;
				v2 = v1;
				v1 = tmp2;
			}
			e.v1(v1);
			e.v2(v2);
			e.c(count);
		}
		else if (param == "dist") {
			tmp >> dist;

			e.dist(dist);
		}
		else if (param == "1direct") {
			tmp >> d1;
			e.d1(d1);
		}
		else if (param == "2direct") {
			tmp >> d2;
			e.d2(d2);

		//	e.debug();
		//	cout << v1 << " " << v2 << " " << count << " " << dist << " " << d1 << " " << d2 << endl;
			//if ((v1 == 20) && (v2 == 1395440)) { 
			//	e.debug();
			//}
			if (count >= minimum_pair) { 
				ep++;

				try {
					add_edge_to_vertex(&verts, v1, e);
					add_edge_to_vertex(&verts, v2, e);
				}
				catch (int x) {
					cout <<"Bad edge at line: " << edgeline << endl;
					e.debug();
				}
			}
			else {
			//	orphan << "Here " << " " << v1 << " " << v2 << " " <<count << " " << minimum_pair << endl;
				orphan << orphan_string;
				orphan_count++;
				ep++;
			}
				
			orphan_string.clear();
			e.init();
		}
	}
	orphan.close();
	

	cerr << "Edges processed: " << ep << endl;
	return verts;
}

void add_edge_to_vertex( map<unsigned int, vertex> *v, unsigned int number,  edge e) { 
	vertex tmp_vertex;
	map <unsigned int, vertex>::iterator v_it;

	v_it = (*v).find(number);
	if (number == 0) { 
		cerr << "Couldn't find vertex: " << number << ". Probably bad sam file. " << endl; 
		e.debug();
		throw 20; 
	}

	(*v_it).second.insert_edge(e);
}

void printUsage() { 
	cerr << "Usage: adjust_graph -g <graph file> -f <contig file> -p <minimum pair count>" << endl;
}






bool vertex::edge_direction(edge e) { 
	//cout << name_private << endl;
	if (name_private == e.v1()) {
		if (e.d1() != 0) { return 1; }
		return 0;
	}
	else if (name_private == e.v2()) {
		if (e.d2() != 0) { return 1; }
		return 0;
	}
	else { cout << "Edge doesn't belong to vertex " << name_private << " "  << e.v1() << " " << e.v2() << endl; throw "Error";}
	
	return 0;
}


bool edge::is_edge(unsigned int x, unsigned int y) { 
	if (y < x) { 
		unsigned int tmp = x;
		x = y;
		y = tmp;
	}

	if ((x == vertex1) && (y == vertex2)) { return 1; }
	return 0;
}

bool edge::is_edge(edge e) { 
	if ((vertex1 == e.v1()) && (vertex2 == e.v2())) { return 1; }
	return 0;
}

bool vertex::is_vertex(unsigned int x) { 
	if (x == name_private) {
		return 1;
	}
	return 0;
}

