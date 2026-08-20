///
///  \file   MMHit.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///

#ifndef MMHit_H
#define MMHit_H



class MMHit : public TObject{
private:
    float x;
    float y;
    float q;
    float yer;
    double fdtime;
    double chnl;
    long m_clusID=-1;

    
public:
    //    void reset();
    
    void SetStrip(float v)    {x=v;}
    void sety(float v)        {y=v;}
    void SetCharge(float v)    {q=v;}
    void setyer(float v)    {yer=v;}
    void SetTime(double t)  {fdtime=t;};
//    void SetChannel(double ch)  const {chnl=ch;};

    void SetClusID(long clusid);

    
    float Strip()    {return (x);}
    float gety()    {return y;}
    float Charge()    {return q;}
    float getyer()    {return yer;}
    float Time()    {return fdtime;}
    long GetClusID() ;


    double Channel()   {return (double)(Strip()); };
    
    void show() {std::cout<<" x = "<<x<<" y = "<<y <<" q = "<<q<<endl;}
    void show_x(){ std::cout<<"x = "<<x<<endl;}
    void show_y(){ std::cout<<"y = "<<y<<endl;}
    void show_q(){ std::cout<<"q = "<<q<<endl;}
    void show_yer(){ std::cout<<"yer = "<<yer<<endl;}
    
    ClassDef(MMHit, 1)
    
};


 void MMHit::SetClusID(long clusid){
  m_clusID = clusid;
}

 long MMHit::GetClusID()  {
  return m_clusID;
}



#endif
