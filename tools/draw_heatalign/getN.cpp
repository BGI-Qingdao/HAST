#include <iostream>
#include <cctype>
#include <string>

void printNZone(const std::string & name , const std::string & seq ) {
    if(name == "" || seq == "" ){
        return ;
    }
    int prev_n = -1 ;
    int curr_n = -1 ;
    for( int i = 0 ; i<(int)seq.size() ; i++ ){
        if( seq.at(i) == 'N' || seq.at(i) == 'n' ) {
            curr_n = i ;
            if(prev_n == -1 ) {
                prev_n = i ;
            }
        } else {
            if( prev_n != -1 && curr_n != -1  ) {
                std::cout<<name<<'\t'<<prev_n+1<<'\t'<<curr_n+1<<'\n';
                prev_n=-1;curr_n=-1;
            }
        }
    }
}

int main(int argc , char ** /**argv**/){
    if( argc != 1 ) {
        std::cerr<<"getN <xxx.fa >xxx.nzone.txt"<<std::endl;
        return 1;
    }
    std::string line;
    std::string name;
    std::string seq;
    while(!std::getline(std::cin,line).eof()){
        if( line.empty() ) continue ;
        if( line[0] == '>' ) {
            printNZone(name,seq);
            name="";
            seq="";
            for(int i = 1 ; i<(int)line.size() ; i++ ){
                if( std::isblank(line[i]) ) break;
                name += line[i];
            }
        }
        else {
            seq += line ;
        }
    }
    printNZone(name,seq);
    return 0;
}
