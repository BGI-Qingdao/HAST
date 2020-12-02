#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <cassert>

/*******************************************************************************
 *
 * Global variables and magic numbers
 *
 *******************************************************************************/

const int total_colors = 10;
const std::string colors1[total_colors] = { "lime","aqua","darkcyan","lawngreen","navy","olive","deepskyblue","indigo","dodgerblue","greenyellow"};
const std::string colors2[total_colors] = { "darkorange","gold","deeppink","crimson","sandybrown","purple","firebrick", "peru","violet","salmon"};

const std::string heatcolors[11] = {
                              "rgba(253,254,191,0.90)" 
                            , "rgba(249,226,123,0.90)"
                            , "rgba(252,191,84, 0.90)"
                            , "rgba(246,159,95, 0.90)"
                            , "rgba(231,133,117,0.90)"
                            , "rgba(207,115,136,0.90)"
                            , "rgba(180,103,149,0.90)"
                            , "rgba(151,93 ,154,0.90)"
                            , "rgba(122,83 ,149,0.90)"
                            , "rgba(92, 85 ,117,0.90)"
                            , "rgba(77, 77 ,79 ,0.90)"
};

const int scale_len = 5000000;
const int scale_step = 200000;
const int scale_lable_step=1000000;
const float min_idy = 0.89;
const std::string ref_name("GRCH38 MHC");
/*******************************************************************************
 *
 * Structure and functions 
 *
 *******************************************************************************/

void PrintUsage(){
    std::cerr<<"Usage :"<<std::endl;
    std::cerr<<"    draw_heatalign ref_len -i H1.align.txt [-i H2.align.txt] ... [-g genes.txt] "<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"Input description   :"<<std::endl;
    std::cerr<<"    xxx.align.txt"<<std::endl;
    std::cerr<<"        filename must follow xxx.align.txt pattern and xxx will print in graph as query name. "<<std::endl;
    std::cerr<<"        content folllow below format :"<<std::endl;
    std::cerr<<"        REF\tREF_START\tREF_END\tQUERY\tQUERY_START\tQUERY_END\t+/-/N\tDIFF"<<std::endl;
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

struct SVG_Align {
    static int graph_width;
    static int graph_height;
    static int ref_len ;
    static int align_num;
    static float scale ;

    static void PrintHeader(){
        std::cout<<R"(<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg height=")"<<graph_height<<R"(" width=")"<<graph_width<<R"(" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">)"<<'\n';
    }

    static void PrintFooter(){
        std::cout<<R"(</svg>)"<<'\n';
    }

    static void Init( int alignnum , int rl){
        align_num=alignnum;
        graph_width = 1200;
        graph_height =  ( int((align_num-1)/2) +1 ) * 120 + 100 ;
        /******************************************************
         *   a pair of alignment --> 100 pixel
         *   upper margin and heatcolor --> 50 pixel
         *   bottom margin and scale    --> 50 pixel
         ******************************************************/
        ref_len = rl;
        scale = float(800) / float(ref_len);
    }

    static float x_pos(int pos) {
        // leave 20 pexel as left margin
        return 50+ float(pos) * scale ;
    }

    static float y_in_ref( int align_index ) {
        return  (((align_index-1)/2) +1)*120;
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

    static float y_scale(){
        return y_in_ref(align_num) + 60;
    }

    //static void PrintTitle(const std::string & title){
    //    std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="165" y="14">)"<<title<<R"(</text>)"<<'\n';
    //}

    static void PrintRefName(const std::string & name, int align_index){
        int y = y_in_ref(align_index) - 6;
        std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="70" y=")"<<y<<R"(">)"<<name<<R"(</text>)"<<'\n';
    }
    static void PrintQueryName(const std::string & name, int align_index){
        int y = y_in_query(align_index);
        if(align_index %2 == 1 )
            y += 15;
        else 
            y -= 6 ;
        std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="70" y=")"<<y<<R"(">)"<<name<<R"(</text>)"<<'\n';
    }

    static void PrintScale(){
        int y = y_scale();
        std::cout<<R"(<line fill="black" stroke="black" stroke-width="1" x1="30" x2="870" y1=")"<<y<<R"(" y2=")"<<y<<R"(" />)"<<'\n';
        for(int pos =0 ; pos <=scale_len; pos += scale_step){
            int x=x_pos(pos);
            int y1 ;
            if( pos % scale_lable_step == 0 )
                y1 = y+5;
            else
                y1 = y+3;

            std::cout<<R"(<line fill="black" stroke="black" stroke-width="1" x1=")"<<x<<R"(" x2=")"<<x<<R"(" y1=")"<<y<<R"(" y2=")"<<y1<<R"(" />)"<<'\n';
            if( pos % scale_lable_step == 0 ) {
                int xx = pos / scale_lable_step;
                std::cout<<R"(<text font-family="Arial" font-size="0.7em" x=")"<<x-10<<R"(" y=")"<<y+15<<R"(">)"<<xx<<R"( Mb </text>)"<<'\n';
            }
        }
    }

    static void PrintGeneText(const std::map<int,std::string> &genes){
        int index = 0;
        int y = y_in_ref(align_num) ;
        for( const auto & pair : genes ){
            index ++ ;
            int x = x_pos(pair.first);
            if( pair.second.size() <3 ){
                if( index %2 ==1 ) {
                    int y1 =y+13 ;
                    std::cout<<R"(<text font-family="Arial" font-size="0.5em" x=")"<<x<<R"(" y=")"<<y1<<R"(" fill="black" >)"<<pair.second<<R"(</text>)"<<'\n';
                } else {
                    int y1 =y-5;
                    std::cout<<R"(<text font-family="Arial" font-size="0.5em" x=")"<<x<<R"(" y=")"<<y1<<R"(" fill="black" >)"<<pair.second<<R"(</text>)"<<'\n';
                }
            } else {
                if( index %2 ==1 ) {
                    int y1 =y+8 ;
                    std::cout<<R"(<text font-family="Arial" font-size="0.5em" x=")"<<x<<R"(" y=")"<<y1<<R"(" fill="black" transform="rotate(60,)"<<x<<','<<y1<<")\">"<<pair.second<<R"(</text>)"<<'\n';
                } else {
                    int y1 =y-5;
                    std::cout<<R"(<text font-family="Arial" font-size="0.5em" x=")"<<x<<R"(" y=")"<<y1<<R"(" fill="black" transform="rotate(-60,)"<<x<<','<<y1<<")\">"<<pair.second<<R"(</text>)"<<'\n';
                }
            }
        }
    }

    static void PrintGenePoint(const std::map<int,std::string> &genes, int align_index){
        for( const auto & pair : genes ){
            PrintPointInRef(pair.first , align_index);
        }
    }

    static void PrintRefLine( int align_index ){
         int y=y_in_ref(align_index);
         std::cout<<"<line fill=\"rgb(112,173,71)\" stroke=\"rgb(112,173,71)\""<<R"( stroke-width="3" x1="50" x2="850" y1=")"<<y<<R"(" y2=")"<<y<<R"(" />)"<<'\n';
    }

    static void PrintBorder(){
        std::cout<<R"(<rect width=")"<<graph_width<<R"(" height=")"<<graph_height<<'"';
        std::cout<<R"( style="fill:rgb(255,255,255);stroke-width:1;stroke:rgb(0,0,0)"<<")\"/>\n";
    }

    static std::string  query_color( int align_index ,int /*color_index*/ ) {
        if( align_index %2 == 1 ) return "rgb(237,125,49)"; // colors1[(align_index/2 + color_index)%total_colors];
        return "rgb(91,155,213)";  // colors2[(align_index/2 + color_index)%total_colors];
    }

    static void PrintPointInRef(int pos , int align_index ){
        int x=x_pos(pos);
        int y=y_in_ref(align_index);
        std::cout<<R"(<circle cx=")"<<x<<R"(" cy=")"<<y<<R"(" r="1" stroke="black" stroke-width="1" fill="black" />)"<<'\n';
    };

    static void PrintPointInQuery(int pos , int align_index){
        int x=x_pos(pos);
        int y=y_in_query(align_index);
        std::cout<<R"(<circle cx=")"<<x<<R"(" cy=")"<<y<<R"(" r="2" stroke="black" stroke-width="2" fill="black" />)"<<'\n';
    };

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
        // TODO :
        if( idy < 0.89 ) color = heatcolors[10];
        if( idy == 1 )  color = heatcolors[0];
        else color = heatcolors[99-int(idy*100)];
        //std::string color = "rgv"
        std::cout<<R"(<polygon points=")"<<x_r_start<<','<<y_ref<<' ';
        std::cout<<x_r_end<<','<<y_ref<<' ';
        std::cout<<x_q_end<<','<<y_query<<' ';
        std::cout<<x_q_start<<','<<y_query<<"\" ";
        std::cout<<R"(style="fill:)"<<color<<R"(;stroke:)"<<color<<R"(;stroke-width:1;" />)"<<'\n';
    }

    //
    static void PrintHeatBar(){
        int y = 15;
        for( int i = 0 ; i<=10 ; i++ ) {
            int x = 100+i*15;
            std::string color=heatcolors[i];
            std::cout<<R"(<rect width=")"<<15<<R"(" height=")"<<15<<R"(" x=")"<<x<<R"(" y=")"<<y<<R"(" )";
            std::cout<<R"(style="fill:)"<<color<<R"(;stroke:)"<<color<<R"(;stroke-width:1;" />)"<<'\n';
        }
        std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="100" y=")"<<45<<R"(">)"<<"0%"<<R"(</text>)"<<'\n';
        std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="250" y=")"<<45<<R"(">)"<<"10%"<<R"(</text>)"<<'\n';
        std::cout<<R"(<text font-family="Arial" font-size="0.7em" x="275" y=")"<<25<<R"(">)"<<"Est.difference"<<R"(</text>)"<<'\n';
    }
};


int SVG_Align::graph_width;
int SVG_Align::graph_height;
int SVG_Align::ref_len ;
int SVG_Align::align_num;
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
            if(o == '+' ) {
                orient = true ;
            } 
            else {
                orient = false ;
                if( query_start < query_end )
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
                 ,ref_pos_max(-1)
    {
    }

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
            SVG_Align::PrintMapRect(block.ref_start,block.ref_end,
                    pos_in_line(block.query_start) ,
                    pos_in_line(block.query_end) ,
                    align_index,
                    block.idy);
        }
        //SVG_Align::PrintPointInQuery(pos_in_line(query_pos_max),align_index);
        //SVG_Align::PrintPointInQuery(pos_in_line(query_pos_min),align_index);
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
    std::vector<bool> n_flags;
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
        int low_idy = 0;
        while( ! std::getline(ifs,line).eof() ){
            AlignBlock block ;
            block.InitFromString(line);
            if(block.idy<min_idy || std::abs(block.ref_end) -std::abs(block.ref_start) < 2000) { low_idy++ ; continue ;}
            if( curr_query_name == ""  || curr_query_name != block.query_name){
                curr_query_name = block.query_name;
                FlushLastSeq();
                seqs.push_back(QuerySeq());
                seqs.at(seqs.size()-1).seq_name = curr_query_name;
            }
            assert(seqs.size()>0);
            seqs.at(seqs.size()-1).AddBlock(block);
        }
        FlushLastSeq();
        for( auto & seq : seqs )
            seq.DetectOrient();
        std::cerr<<"filter "<<low_idy<<" low idy maps by min_idy="<<min_idy<<std::endl;
        std::cerr<<"loading data end with "<<seqs.size()<<" query sequence(s)."<<std::endl;
        ifs.close();
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

void LoadGenes(const std::string & fn , std::map<int, std::string> & gs ) {
    std::ifstream ifs(fn);
    if( ! ifs.is_open() ){
        std::cerr<<"failed to open file :"<<fn<<std::endl;
        std::cerr<<"exit ... "<<std::endl;
        exit(1);
    }
    std::cerr<<"loading data from "<<fn<<std::endl;
    std::string line;
    while( ! std::getline(ifs,line).eof() ){
        int pos ; std::string name;
        std::istringstream ist(line);
        ist>>pos>>name;
        gs[pos]=name;
    }
    ifs.close();
    std::cerr<<"load "<<gs.size()<<" genes from "<<fn<<std::endl;
}

/*******************************************************************************
 *
 * main function
 *
 *******************************************************************************/

int main( int argc , char ** argv ) {
    if( argc <4 || argc%2!=0 ){
        PrintUsage();
        return 1;
    }
    // caches
    std::string gene_file;
    std::vector<std::pair<std::string, std::string> > aligns;
    std::vector<Query> querys;
    std::map<int,std::string> genes;
    // parse parameters
    int draw_len = std::atoi(argv[1]);
    for( int i = 2; i < argc-1 ; i+=2 ) {
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
    // load files into caches
    int align_num= 0;
    for( const auto & pair : aligns ){
        align_num ++ ;
        Query q ;
        q.query_name = pair.first ;
        q.align_index = align_num;
        q.LoadAlignFile(pair.second);
        querys.push_back(q);
    }
    if( gene_file != "" ) {
        LoadGenes(gene_file,genes);
    }
    //print SVG:
    SVG_Align::Init(align_num, draw_len);
    SVG_Align::PrintHeader();
    SVG_Align::PrintBorder();
    SVG_Align::PrintHeatBar();

    for( const auto & query : querys ) {
        if(query.align_index %2 == 1 ){
            SVG_Align::PrintRefLine(query.align_index);
            if( gene_file != "" ) {
                SVG_Align::PrintGenePoint(genes,query.align_index);
            }
        }
        query.PrintQueryLine();
        query.PrintAlignRect();
        if(query.align_index %2 == 1 ){
            SVG_Align::PrintRefName(ref_name,query.align_index);
        }
        SVG_Align::PrintQueryName(query.query_name ,query.align_index);
    }
    SVG_Align::PrintScale();
    if( gene_file != "" ) {
        SVG_Align::PrintGeneText(genes);
    }
    SVG_Align::PrintFooter();
    // done
    return 0;
}
