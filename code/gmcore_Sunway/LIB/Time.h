#ifndef TIME_H_INCLUDED
#define TIME_H_INCLUDED 1

typedef struct{
  int oldt;
  int newt;
  int timestep;
}TimeType; 

TimeType Timeinfo;

void TimeInit();
void Time_Advance();
#endif