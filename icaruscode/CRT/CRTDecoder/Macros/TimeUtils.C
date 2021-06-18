#define TimeUtils_C

struct Date{
        int year, month, day;
};

const int monthDays[12] = {31, 28, 31, 30, 31, 30,
                           31, 31, 30, 31, 30, 31};

int MonthNameToNum(string name){
        if(name=="jan") return 1;
        if(name=="feb") return 2;
        if(name=="mar") return 3;
        if(name=="apr") return 4;
        if(name=="may") return 5;
        if(name=="jun") return 6;
        if(name=="jul") return 7;
        if(name=="aug") return 8;
        if(name=="sep") return 9;
        if(name=="oct") return 10;
        if(name=="nov") return 11;
        if(name=="dec") return 12;
        cout << "month not identified!" << endl;
        return -1;
}

Date PathToDate(string path){
        string ymd=path.substr(path.find_last_of("/")+1);
        ymd = ymd.substr(0,ymd.find("_"));
        if(ymd.size()!=9)
                cout << "WARNING: file name format error (" << ymd << ")" << endl;

        Date date;
        string tmp = ymd.substr(0,4);
        date.year = atoi(tmp.c_str());

        tmp = ymd.substr(4,3);
        date.month = MonthNameToNum(tmp);

        tmp = ymd.substr(7);
        date.day = atoi(tmp.c_str());

        return date;
}

// This function counts number of leap years before the 
// given date 
int countLeapYears(Date d)
{
    int years = d.year;

    // Check if the current year needs to be considered 
    // for the count of leap years or not 
    if (d.month <= 2)
        years--;

    // An year is a leap year if it is a multiple of 4, 
    // multiple of 400 and not a multiple of 100. 
    return years / 4 - years / 100 + years / 400;
}

int DaysPassed(Date dt1, Date dt2){
   // COUNT TOTAL NUMBER OF DAYS BEFORE FIRST DATE 'dt1' 

    // initialize count using years and day 
    long int n1 = dt1.year*365 + dt1.day;

    // Add days for months in given date 
    for (int i=0; i<dt1.month - 1; i++)
        n1 += monthDays[i];

    // Since every leap year is of 366 days, 
    // Add a day for every leap year 
    n1 += countLeapYears(dt1);

    // SIMILARLY, COUNT TOTAL NUMBER OF DAYS BEFORE 'dt2' 

    long int n2 = dt2.year*365 + dt2.day;
    for (int i=0; i<dt2.month - 1; i++)
        n2 += monthDays[i];
    n2 += countLeapYears(dt2);

    // return difference between two counts 
    return (n2 - n1);
}
