#ifndef KMER_KMER_H
#define KMER_KMER_H
//#include <stdint.h>
#include <vector>
#include <string>
#include <cassert>

typedef unsigned long long  ubyte8;

struct BaseStr{
    static char base2int(char base) { return  (char)(((base)&0x06)>>1)  ; }   //base ACTG => int 0123
    static char int2base(int seq) { return  "ACTG"[seq] ;}
    static char int_comp(char seq) { return  (char)(seq^0x02) ; }
    static std::string BaseStr2Str(const std::vector<char>& base) {
        int len = base.size();
        std::string ret ;
        if ( len < 1 )
        {
            return ret ;
        }
        ret.resize(len);
        for( int i = 0 ; i < len ; i++ )
            ret[i] = int2base(base.at(i));
        return ret ;
    }
    static std::vector<char> str2BaseStr(const std::string & read ){
        int len = read.size();
        std::vector<char> ret ;
        if ( len < 1 )
        {
            return ret ;
        }
        ret.resize(len);
        for( int i = 0 ; i < len ; i++ )
            ret[i] = base2int(read.at(i));
        return ret ;
    }
    static std::vector<char> reverseComplementSeq ( const std::vector<char> & str  )
    {
        int len = str.size();
        std::vector<char> ret ;
        if ( len < 1 )
        {
            return ret ;
        }
        ret.resize(len);
        for ( int i = len - 1 , index = 0; i >= 0; i-- , index ++ )
        {
            ret[index] = int_comp (str.at(i));
        }
        return ret;
    }
};
struct Kmer
{
    static Kmer WORDFILTER;
    static int overlap;
    static void InitFilter(int toverlap){
        WORDFILTER = createFilter(toverlap);
        overlap = toverlap;
    }

    unsigned long long high, low;
    bool operator < ( const Kmer & kmer1) const
    {
        if ( high < kmer1.high )
        {
            return 1;
        }
        else if ( high == kmer1.high )
        {
            if ( low < kmer1.low )
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 0;
        }
    }
    bool operator == ( const Kmer & kmer1) const
    {
        return high == kmer1.high && low == kmer1.low ;
    }
    void And(const Kmer & kmer2){
        high &= kmer2.high;
        low &= kmer2.low;
    }
    void LeftBitMoveBy2() 
    {
        ubyte8 temp = low >> 62;
        high <<= 2;
        high |= temp;
        low <<= 2;
    }
    void RightBitMoveBy2()
    {
      ubyte8 temp = ( high & 0x3 ) << 62;
      high >>= 2;
      low >>= 2;
      low |= temp;
    }

    void nextKmer (char ch )
    {
        LeftBitMoveBy2();
        And(WORDFILTER );
        low |= ch;
    }

    void prevKmer ( char ch )
    {
        RightBitMoveBy2();
        if ( 2 * ( overlap - 1 ) < 64 )
        {
            low |= ( ( ( ubyte8 ) ch ) << 2 * ( overlap - 1 ) );
        }
        else
        {
            high |= ( ( ubyte8 ) ch ) << ( 2 * ( overlap - 1 ) - 64 );
        }
    }

    static Kmer createFilter ( int overlaplen )
    {
        Kmer word;
        word.high = word.low = 0;

        if ( 2 * overlaplen < 64 )
        {
            word.low = ( ( ( ubyte8) 1 ) << ( 2 * overlaplen ) ) - 1;
        }
        else
        {
            word.low = ~word.low;

            if ( 2 * overlaplen > 64 )
            {
                word.high = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 64 ) ) - 1;
            }
        }
        return word;
    }
    void Init() {
        high = 0 ; low = 0;
    }
    // make sure str is base2int str , not  AGCT str
    static Kmer str2Kmer(const std::vector<char> & str){
        assert(str.size() == overlap);
        Kmer word; word.Init();
        for (int index = 0; index < overlap; index++ )
        {
            word.LeftBitMoveBy2();
            word.low |= str.at(index);
        }
        Kmer bal_word = word.reverseComplement();
        if( word < bal_word ) 
            return word ;
        else
            return bal_word;
    }

    // make sure read is base2int str , not  AGCT str
    static std::vector<Kmer> chopRead2Kmer( const std::vector<char> & read ){
        int rlen = read.size() ;
        assert(rlen >= overlap);
        std::vector<Kmer>  ret;
        Kmer word; word.Init();
        std::vector<char> bal_read = BaseStr::reverseComplementSeq(read);
        for (int index = 0; index < overlap; index++ )
        {
            word.LeftBitMoveBy2();
            word.low |= read.at(index);
        }
        Kmer bal_word = word.reverseComplement();
        if( word < bal_word ) 
            ret.push_back(word) ;
        else
            ret.push_back(bal_word);
        for( int index = 1 ; index <= rlen - overlap ; index ++ ){
            word.nextKmer (read.at(index - 1 + overlap));
            bal_word.prevKmer( bal_read[rlen - index - overlap] );
            if( word < bal_word ) 
                ret.push_back(word) ;
            else
                ret.push_back(bal_word);
        }
        return ret;
    }

    static Kmer fastReverseComp ( const Kmer &base , char seq_size )
    {
        Kmer seq = base ;
        seq.low ^= 0xAAAAAAAAAAAAAAAALLU;
        seq.low = ( ( seq.low & 0x3333333333333333LLU ) << 2 ) | ( ( seq.low & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
        seq.low = ( ( seq.low & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.low & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
        seq.low = ( ( seq.low & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.low & 0xFF00FF00FF00FF00LLU ) >> 8 );
        seq.low = ( ( seq.low & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.low & 0xFFFF0000FFFF0000LLU ) >> 16 );
        seq.low = ( ( seq.low & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.low & 0xFFFFFFFF00000000LLU ) >> 32 );

        if ( seq_size < 32 )
        {
            seq.low >>= ( 64 - ( seq_size << 1 ) );
            return seq;
        }

        seq.high ^= 0xAAAAAAAAAAAAAAAALLU;
        seq.high = ( ( seq.high & 0x3333333333333333LLU ) << 2 ) | ( ( seq.high & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
        seq.high = ( ( seq.high & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.high & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
        seq.high = ( ( seq.high & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.high & 0xFF00FF00FF00FF00LLU ) >> 8 );
        seq.high = ( ( seq.high & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.high & 0xFFFF0000FFFF0000LLU ) >> 16 );
        seq.high = ( ( seq.high & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.high & 0xFFFFFFFF00000000LLU ) >> 32 );
        ubyte8 temp = seq.high;
        seq.high = seq.low;
        seq.low = temp;
        seq.RightBitMove ( 128 - ( seq_size << 1 ) );
        return seq;
    }

    void RightBitMove (int dis )
    {
        if ( dis < 64 )
        {
            ubyte8 mask = ( ( ( ubyte8 ) 1 ) << dis ) - 1;
            ubyte8 temp = ( high & mask ) << ( 64 - dis );
            high >>= dis;
            low >>= dis;
            low |= temp;
        }
        high >>= ( dis - 64 );
        low = high;
        high = 0;
    }

    Kmer reverseComplement ()
    {
        return fastReverseComp (*this, overlap );
    }
    static std::vector<char> ToBaseStr(const Kmer & k){
        Kmer tmp = k;
        std::vector<char> kmer;
        kmer.resize(overlap);
        for( int i = 0 ; i <overlap ; i++ ) {
            char next = tmp.low & 0x3 ;
            kmer[overlap-1-i]=next;
            tmp.RightBitMoveBy2();
        }
        return kmer;
    }
};

namespace std{
	template<>
		struct hash<Kmer>
		{
			size_t operator()(const Kmer& x) const
			{
				return x.low;
			}
		};
}
#endif
