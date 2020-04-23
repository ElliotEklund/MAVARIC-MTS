#include <iostream>

//#define TEST

int main(){

    #ifdef TEST
    std::cout << "Hi" << std::endl;
    #endif
    return 0;
}
