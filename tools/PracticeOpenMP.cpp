#include<stdio.h>
#include<iostream>
#include<omp.h>
using namespace std; 

int main()
{   
    #pragma omp parallel
    for (int i = 0; i < 2; i++)
    {
        cout<< i<< "\n";
    }
}

