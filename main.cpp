/*
C++ code to generate the longest string of length N that can have the assembly index equal to N-1
as a function of the radix of the initial assembly pool

(c) Szymon Lukaszyk
licensed under MIT License
email: szymon@patent.pl
History
v1: 24.09.2024 1st working version

Based on
https://www.mdpi.com/2227-7390/12/10/1600
https://www.preprints.org/manuscript/202409.1581
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include <string>
#include <bitset>
#include <cmath>

using namespace std;

// Convert number 0,1,...,9,10,11,... to string 0,1,...,9,a,b,...
string intToChar(const int &num) {
    if (num >= 10 && num <= 35) {
        // Convert 10 -> 'a', 11 -> 'b', ..., 35 -> 'z'
        char result = 'a' + (num - 10);
        return string(1, result); // Return as a string
    }
    else {
        return to_string(num);
    }
}

// Count occurrences of string subStr within the string& str
int countOccurrences(const string& str, const string& subStr) {
    if (subStr.empty()) {
        return 0;  // Avoid infinite loop if the substring is empty
    }
    int count = 0;
    size_t pos = str.find(subStr);
    // Keep finding the substring until the end of the string
    while (pos != string::npos) {
        count++;
        pos = str.find(subStr, pos + subStr.length());  // Move position ahead by the length of the substring
    }
    return count;
}

// Return the number of repeated doublet dbl within the string& str
int areRepeatedDoublets(const string &str, const string &dbl){
    int count = 0;
    for (size_t i = 0; i < str.size()-1; i++){
        // Check if the current pair matches the doublet
        if( str[i] == dbl[0] && str[i+1] == dbl[1] ){
            count++;
            // Account for triplets 000, etc.
            if( dbl[0] == dbl[1] && str[i+2] == dbl[1]){
                i++;
            }
        }
    }
    return count;
}

// Check if the string str of length N is the longest that can have the assembly index equal to N-1 for b
// Warning: Maximal number of triplets is not checked!
bool checkTheNmax(const string &str, const size_t &b){
    bool retval = true;
    // Generate matrix of doublets
    for (size_t r=0; r<b; r++ ){
        for (size_t c=0; c<b; c++ ){
            string dbl = intToChar(r) + intToChar(c);
            int ret = areRepeatedDoublets(str, dbl);
            if (ret==0) {
                cout << "ERR Doublet " << dbl << " does not appear in the string!" << endl;
                retval = false;
            }
            else if (ret>1) {
                cout << "ERR Doublet " << dbl << " repeats in the string " << ret << " times!" << endl;
                retval = false;
            }
//            else {
//                cout << "OK! Doublet " << dbl << " appears and does not repeat in the string!" << endl;
//            }
        }
    }
    for(size_t k=0; k<b; k++){
        cout << "b=" << intToChar(k) << " occurs " << countOccurrences(str, intToChar(k)) << " times" << endl;
    }

    return retval;
}

// Generate the string str of length N is the longest that can have the assembly index equal to N-1 for b (WB version)
string generateTheNmaxWB(const size_t &b){
    if( b < 1)
        return "ERR b < 1";
    string str = "";
    // base triplets
    for(size_t c = 0; c < b; c++) {
        for(size_t k = 0; k < 3; k++) {
            str = str + intToChar(c);
        }
    }
    if( b == 1)
        return str;
    if( b == 2)
        return str + "0";
    if( b == 3)
        return str + "1020";
    if( b == 4)
        return str + "102132030";
    if( b == 5)
        return str + "1021324303140420";

    for(size_t c = 0; c <= b-2; c++) { //1st pass, 0-4 for b=6, etc.
        str = str + intToChar(c+1) + intToChar(c);
    }
    for(size_t c = 0; c <= b-4; c++) { //2nd pass, 0-2 for b=6, etc.
        str = str + intToChar(c) + intToChar(c+3);
    }
    for(size_t c = 0; c <= b-5; c++) { //3rd pass, 0-1 for b=6, etc.
        str = str + intToChar(c) + intToChar(c+4);
    }
    str = str + "20" + intToChar(b-1) + intToChar(b-3); //4th pass "2053" for b=6, etc.

    size_t npass = b-5; // n. of passes
    if( npass%2 != 0)   // if npass is odd
        npass--;        // round it down to integer
    // https://oeis.org/A263883
    for(size_t pn = 1; pn <= npass; pn++) {
        if( pn%2 != 0){ // odd pass, upper triangle
            size_t start_r = pn+4; // 5,7,9...
            if(start_r+1 < b){
                size_t nel = b-start_r; // n. of elements
                for(size_t c = 0; c <nel; c++) {
                    str = str + intToChar(c) + intToChar(start_r+c);
                }
            }
        }
        else{           // even pass, bottom triangle
            size_t start_r = (pn+4)/2; // 3,4,5,6,7,...
            str = str + intToChar(start_r) + intToChar(0); // 30, 40, 50,...
            size_t next_r = start_r*2; // 6, 8, 10, 12,...
            if( next_r+1 < b){
                size_t nel = b-next_r;
                for(size_t c = 1; c <= nel; c++) {
                    int rel = c;
                    if(c == nel){ //last entry
                        rel = c+start_r-1;
                    }
                    str = str + intToChar(next_r+c-1) + intToChar(rel);
                }
            }
        }
    }
    if( b%2 == 0)   // b is even
        str = str + intToChar(0);
    return str;
}

// Generate the string str of length N is the longest that can have the assembly index equal to N-1 for b (PM version)
string generateTheNmaxPM(const size_t &b){
    if( b < 1)
        return "ERR b < 1";
    string str = "";
    // base triplets
    for(size_t c = 0; c < b; c++) {
        for(size_t k = 0; k < 3; k++) {
            str = str + intToChar(c);
        }
    }
    if( b == 1)
        return str;
    if( b == 2)
        return str + "0";
    if( b == 3)
        return str + "0210";
    if( b == 4)
        return str + "031021320";
    if( b == 5)
        return str + "0410213243031420";
    if( b == 6)
        return str + "0510213243540314253041520";
    if( b == 7)
        return str + "061021324354650314253640516204152630";
    // 1. top right corner
    str = str + "0" + intToChar(b-1);
    // 2. first subdiagonal
    for(size_t c = 1; c <= b-1; c++) {
        str = str + intToChar(c)  + intToChar(c-1);
    }
    // 3. first superdiagonal
    for(size_t c = 0; c <= b-4; c++) {
        str = str + intToChar(c)  + intToChar(c+3);
    }
    // 4. second subdiagonal
    for(size_t c = 3; c <= b-1; c++) {
        str = str + intToChar(c)  + intToChar(c-3);
    }
    // 5. second superdiagonal
    for(size_t c = 0; c <= b-6; c++) {
        str = str + intToChar(c)  + intToChar(c+5);
    }
    size_t npass = b-7; // n. of passes
    for(size_t pn = 6; pn < npass+6; pn++) {
        if( pn%2 != 0){ // 7,9,11,... odd pass, upper triangle
            if(pn+1 < b){
                size_t nel = b-pn+(b%2); // n. of elements
                for(size_t c = 0; c <nel; c++) {
                    //06,17,28,...; 08,19,... for odd b
                    //07,18,29,...; 09,1a,... for even b
                    str = str + intToChar(c) + intToChar(c+pn-(b%2));
                }
            }
        }
        else{ // 6,8,10,... even pass, lower triangle
            if(pn+1 < b){
                size_t nel = b-pn+1-(b%2); // n. of elements
                for(size_t c = 0; c <nel; c++) {
                    str = str + intToChar(c+pn-1+(b%2)) + intToChar(c);
                }
            }
        }
    }
    if( b%2 == 0 ){ // b is even
        str = str + "0";
    }
    else{ // b is odd
        if(b==13)
            str = str + "70";
        else
            str = str + intToChar((b-1)/2) + "0";
    }
    return str;
}

int main() {
    string str;
    int b;

    ofstream outFile("output.txt");
    if (outFile.is_open()) {
        for(b=1; b<36; b++){
        //for(b=13; b<14; b++){
            str = generateTheNmaxWB(b);
            //str = generateTheNmaxPM(b);
            cout << "b = " << b << " Length: " << str.size()  << endl;
            cout << "String: " << str << endl;
            outFile << "b = " << b << " Length: " << str.size()  << endl;
            outFile << "String: " << str << endl;
            checkTheNmax(str, b);
        }
    }
    outFile.close();
    return 0;
}
