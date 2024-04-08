#include<stdio.h>

#include "Time.h"

void TimeInit(){
    Timeinfo.oldt = 0;
    Timeinfo.newt = 1;
    Timeinfo.timestep = 0;
}

void Time_Advance(){
    int tmp;

    Timeinfo.oldt ^= 1;
    Timeinfo.newt ^= 1;
    Timeinfo.timestep += 1;
}