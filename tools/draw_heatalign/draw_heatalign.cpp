#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <cassert>

void PrintUsage(){
    std::cerr<<"Usage :"<<std::endl;
    std::cerr<<"    draw_heatalign graph_name ref_len ref_shift -i H1.align.txt [-i H2.align.txt] ... [-g genes.txt] "<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"Input description   :"<<std::endl;
    std::cerr<<"    xxx.align.txt"<<std::endl;
    std::cerr<<"        filename must follow xxx.align.txt pattern and xxx will print in graph as query name. "<<std::endl;
    std::cerr<<"        content folllow below format :"<<std::endl;
    std::cerr<<"        REF\tREF_START\tREF_END\tQUERY\tQUERY_START\tQUERY_END\tDIFF"<<std::endl;
    std::cerr<<"        each line represent a alignment block and blocks from 1 query sequence must stay adjacent."<<std::endl;
    std::cerr<<"    genes.txt"<<std::endl;
    std::cerr<<"        filename has no limit but contents need folllow below format :"<<std::endl;
    std::cerr<<"        position\tname"<<std::endl;
    std::cerr<<std::endl;
    return ;
}

bool parse_align_file_name(const std::string & line, std::string & name) {
    if(line.size() < 11 || line.substr(line.size()-10,10) != ".align.txt" ) {
        return false ;
    }
    name=line.substr(0,line.size()-10);
    return true;
}

int total_colors = 4;
std::string colors1[4] = { "lime","aqua","lawngreen","greenyellow" };
std::string colors2[4] = { "orange","orangered","sandybrown","salmon" };

std::string heatcolors[11] = {"rgba(253,254,191,0.75)" 
                            , "rgba(249,226,123,0.75)"
                            , "rgba(252,191,84, 0.75)"
                            , "rgba(246,159,95, 0.75)"
                            , "rgba(231,133,117,0.75)"
                            , "rgba(207,115,136,0.75)"
                            , "rgba(180,103,149,0.75)"
                            , "rgba(151,93 ,154,0.75)"
                            , "rgba(122,83 ,149,0.75)"
                            , "rgba(92, 85 ,117,0.75)"
                            , "rgba(77, 77 ,79,0.75)"
};

struct SVG_Align {
    static int graph_width;
    static int graph_height;
    static int ref_len ;
    static float scale ;

    static void PrintHeader(){
        std::cout<<R"(<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg height=")"<<graph_height<<R"(" width=")"<<graph_width<<R"(" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">)"<<'\n';
    }

    static void PrintFooter(){
        std::cout<<R"(</svg>)"<<'\n';
    }

    static void Init( int align_num , bool gene , int rl){
        graph_width = 900;
        graph_height = 50 + ( int((align_num-1)/2) +1 ) * 100 + 100 +(gene ? 50 : 0);
        /******************************************************
         *   margin --> 10 pixel
         *   title --> 30 pixel
         *   scale --> 9 pixel
         *   a pair of alignment --> 50 pixel
         *      5 pixel query * 2
         *      6 pixel ref
         *      17 pixel aligned rect * 2
         *
         *   gene names --> 50 pixel
         *
         ******************************************************/
        ref_len = rl;
        scale = float(800) / float(ref_len);
    }

    static float x_pos(int pos) {
        // leave 20 pexel as left margin
        return 50+ float(pos) * scale ;
    }

    static float y_in_ref( int align_index ) {
        return 50 + (((align_index-1)/2) +1)*100;
    }

    static float y_in_ref_rect(int align_index ) {
        if( align_index % 2 == 1) 
            return y_in_ref(align_index)-2;
        return y_in_ref(align_index)+2;
    }
    static float y_in_query_rect(int align_index ) {
        if( align_index % 2 == 1) 
            return y_in_ref(align_index)-45;
        return y_in_ref(align_index)+45;
    }

    static float y_in_query(int align_index ) {
        if( align_index % 2 == 1) 
            return y_in_ref(align_index)-47;
        return y_in_ref(align_index)+47;
    }

    static void PrintTitle(const std::string & title){
        std::cout<<R"(<text font-family="TimeNewRoman" font-size="0.7em" x="165" y="14">)"<<title<<R"(</text>)"<<'\n';
    }

    static void PrintScale(int ref_shift){

    }

    static void PrintGeneText(const std::map<int,std::string> &genes, int align_index){

    }
    static void PrintGenePoint(const std::map<int,std::string> &genes, int align_index){

    }

    static void PrintRefLine( int align_index ){
         int y=y_in_ref(align_index);
         std::cout<<R"(<line fill="blue" stroke="blue" stroke-width="3" x1="50" x2="850" y1=")"<<y<<R"(" y2=")"<<y<<R"(" />)"<<'\n';
    }

    static void PrintBorder(){
        std::cout<<R"(<rect width=")"<<graph_width<<R"(" height=")"<<graph_height<<'"';
        std::cout<<R"( style="fill:rgb(255,255,255);stroke-width:1;stroke:rgb(0,0,0)"<<")\"/>\n";
    }

    static const std::string & query_color( int align_index ,int color_index ) {
        if( align_index %2 == 0 ) return colors1[(align_index/2 + color_index)%total_colors];
        return colors2[color_index%total_colors];
    }

    static void PrintQueryLine(int start , int end, int align_index ,int color_index) {
        int x1 = x_pos(start) ;int x2 = x_pos(end) ;
        int y = y_in_query(align_index);
        const std::string & color = query_color(align_index,color_index);
        std::cout<<R"(<line fill=")"<<color<<R"(" stroke=")"<<color<<R"(" stroke-width="3")";
        std::cout<<R"( x1=")"<<x1<<R"(" x2=")"<<x2<<R"(" y1=")"<<y<<R"(" y2=")"<<y<<R"(" />)"<<'\n';
    }

    static void PrintMapRect(int rstart , int rend , int qstart , int qend, int align_index,float idy){
        float x_r_start = x_pos(rstart) ;
        float x_r_end = x_pos(rend);
        float x_q_start = x_pos(qstart);
        float x_q_end = x_pos(qend);
        float y_ref = y_in_ref_rect(align_index);
        float y_query = y_in_query_rect(align_index);
        std::string color;
        if( idy < 0.8 ) color = heatcolors[10];
        else color = heatcolors[int((1.0-idy)*50)];
        //std::string color = "rgv"
        std::cout<<R"(<polygon points=")"<<x_r_start<<','<<y_ref<<' ';
        std::cout<<x_r_end<<','<<y_ref<<' ';
        std::cout<<x_q_end<<','<<y_query<<' ';
        std::cout<<x_q_start<<','<<y_query<<"\" ";
        std::cout<<R"(style="fill:)"<<color<<R"(;stroke:none;stroke-width:1;" />)"<<'\n';
    }

    static void PrintQueryName( const std::string & name , int align_index){

    }
};


int SVG_Align::graph_width;
int SVG_Align::graph_height;
int SVG_Align::ref_len ;
float SVG_Align::scale ;


struct AlignBlock {
    std::string m_line ;
    int ref_start;
    int ref_end ;
    std::string ref_name;
    std::string query_name;
    int query_start;
    int query_end;
    float idy;
    bool orient;
    int maped_len() const {
        return ref_end - ref_start + 1 ;
    }
    void InitFromString(const std::string & line){
        m_line=line;
        int det=0;
        for(char c:line) if(c=='\t')det++;
        if(det<6) {
            std::cerr<<"align info is invalid:"<<std::endl;
            std::cerr<<line<<std::endl;
            std::cerr<<"please use \\t to seperate columns!!!"<<std::endl;
            std::cerr<<"exit ..."<<std::endl;
        }
        std::istringstream ist(line);
        if( det == 6 ){
            ist>>ref_name>>ref_start>>ref_end>>query_name>>query_start>>query_end>>idy;
            if( query_start < query_end ) orient = true;
            else orient = false ;
        }
        if( det >= 7 ) {
            char o;
            ist>>ref_name>>ref_start>>ref_end>>query_name>>query_start>>query_end>>o>>idy;
            assert(query_start < query_end );
            if(o == '+' ) {
                orient = true ;
            } 
            else {
                orient = false ;
                std::swap(query_start,query_end);
            }
        }
    }
};
void update_min( int & min , int curr ) {
    if( min == -1 || min > curr ) min = curr;
}
void update_max( int & max , int curr ) {
    if( max == -1 || max < curr ) max = curr;
}

struct QuerySeq {

    int query_shift;
    int query_pos_min;
    int query_pos_max;
    int ref_pos_min;
    int ref_pos_max;
    std::string seq_name;
    QuerySeq() : query_shift(0)
                 ,query_pos_min(-1)
                 ,query_pos_max(-1)
        ,ref_pos_min(-1)
                 ,ref_pos_max(-1){}

        int seq_len() const {
            return query_pos_max - query_pos_min + 1;
        }
        int line_start() const {
            return query_shift ;
        }
        int line_end() const {
            return query_shift + seq_len();
        }
        // Input pos : pos in query seq. 1 base
        // Output    : pos in draw line. 0 base
        float pos_in_line( int pos ) const {
            if( orient )
                return (pos - query_pos_min + query_shift);
            else 
                return (query_shift + query_pos_max - pos );
        }

        std::vector<AlignBlock> blocks;

        void SetShift(int prev_line_end) {
            assert(blocks.size() > 0);
            for( const auto & b : blocks ) {
                update_min(query_pos_min, b.query_end);
                update_min(query_pos_min, b.query_start);
                update_max(query_pos_max, b.query_start);
                update_max(query_pos_max, b.query_end);

                update_min(ref_pos_min, b.ref_end);
                update_min(ref_pos_min, b.ref_start);
                update_max(ref_pos_max, b.ref_start);
                update_max(ref_pos_max, b.ref_end);
            }
            if( prev_line_end < ref_pos_min ) 
                query_shift= ref_pos_min;
            else 
                query_shift =  prev_line_end;
            std::cerr<<"    --> load query "<<seq_name<<" with "<<blocks.size()<<" blocks in "<<seq_len()<<" bps."<<std::endl;
        }

        void AddBlock(const AlignBlock & b) {
            blocks.push_back(b);
        }

        void PrintAlignRect(int align_index) const {
            for(const auto & block : blocks) {
                if( block.maped_len() < 500 ) continue ;
                SVG_Align::PrintMapRect(block.ref_start,block.ref_end,
                        pos_in_line(block.query_start) ,
                        pos_in_line(block.query_end) ,
                        align_index,
                        block.idy);
            }
        }
        bool orient;
        void DetectOrient() {
            int t = 0 , f = 0 ;
            for(const auto & block : blocks) {
                if(block.orient) t+=block.maped_len();
                else f+=block.maped_len();
            }
            if( t > f ) orient = true ;
            else orient = false ;
        }
    };

struct Query {
    std::string query_name ;
    int align_index;
    std::vector<QuerySeq> seqs;
    void FlushLastSeq(){
        if( seqs.size() == 1) {
            seqs.at(0).SetShift(0);
        } else if ( seqs.size() >1 ) {
            int prev_line_end = seqs.at(seqs.size()-2).line_end();
            seqs.at(seqs.size()-1).SetShift(prev_line_end);
        }
    }
    void  LoadAlignFile(const std::string & filename) {
        std::ifstream ifs(filename);
        if( ! ifs.is_open() ){
            std::cerr<<"failed to open file :"<<filename<<std::endl;
            std::cerr<<"exit ... "<<std::endl;
            exit(1);
        }
        std::cerr<<"loading data from "<<filename<<std::endl;
        std::string line;
        std::string curr_query_name;
        while( ! std::getline(ifs,line).eof() ){
            AlignBlock block ;
            block.InitFromString(line);
            if( curr_query_name == ""  || curr_query_name != block.query_name){
                curr_query_name = block.query_name;
                FlushLastSeq();
                seqs.push_back(QuerySeq());
                seqs.at(seqs.size()-1).seq_name = curr_query_name;
            }
            else {
                assert(seqs.size()>0);
                seqs.at(seqs.size()-1).AddBlock(block);
            }
        }
        FlushLastSeq();
        for( auto & seq : seqs )
            seq.DetectOrient();
        std::cerr<<"loading data end with "<<seqs.size()<<" query sequence(s)."<<std::endl;
    }
    void PrintQueryLine() const {
        int color_index = 0 ;
        for(const auto & seq :seqs ){
            SVG_Align::PrintQueryLine(seq.line_start() , seq.line_end() ,align_index,color_index);
            color_index++;
        }
    }

    void PrintAlignRect() const {
        for(const auto & seq :seqs ){
            seq.PrintAlignRect(align_index);
        }
    }
};

void  LoadGenes(const std::string & f) {

}

std::string gene_file;
std::vector<std::pair<std::string, std::string> > aligns;
std::vector<Query> querys;
std::map<int,std::string> genes;

int main( int argc , char ** argv ) {
    if( argc <6 || argc%2!=0 ){
        PrintUsage();
        return 1;
    }
    std::string graph_name = std::string(argv[1]);
    int draw_len = std::atoi(argv[2]);
    int ref_shift = std::atoi(argv[3]);
    for( int i = 4; i < argc-1 ; i+=2 ) {
        if( std::string(argv[i]) == "-i" ){
            std::string name;
            if( parse_align_file_name( std::string(argv[i+1]) ,name ) ){
                aligns.push_back(std::make_pair(name,std::string(argv[i+1])));
            } 
            else {
                PrintUsage();
                return 1;
            }
        }
        else if ( std::string(argv[i]) == "-g" ){
            if(gene_file != "" ){
                std::cerr<<" ERROR : only support one \"-g gene.txt\" but multi -g detected !"<<std::endl;
                PrintUsage();
                return 1;
            }
            gene_file = std::string(argv[i+1]);
        }
        else {
            PrintUsage();
            return 1;
        }
    }
    int align_index  = 0;
    for( const auto & pair : aligns ){
        align_index ++ ;
        Query q ;
        q.query_name = pair.first ;
        q.align_index = align_index ;
        q.LoadAlignFile(pair.second);
        querys.push_back(q);
    }
    if( gene_file != "" ) {
        LoadGenes(gene_file);
    }
    //Draw :
    SVG_Align::Init(align_index,!gene_file.empty() , draw_len);
    SVG_Align::PrintHeader();
    SVG_Align::PrintBorder();
    SVG_Align::PrintTitle(graph_name);

    for( const auto & query : querys ) {
        if(query.align_index %2 == 1 ){
            SVG_Align::PrintRefLine(query.align_index);
            SVG_Align::PrintGenePoint(genes,query.align_index);
        }
        query.PrintQueryLine();
        query.PrintAlignRect();
    }
    SVG_Align::PrintScale(ref_shift);
    SVG_Align::PrintFooter();
    return 0;
}
