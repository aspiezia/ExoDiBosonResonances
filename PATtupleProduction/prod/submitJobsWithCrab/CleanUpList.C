#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TString.h"
typedef struct{
  TString name;
  int idx;
  int subidx;
  long long size;
} filedata;


bool decision(filedata f1, filedata f2){
  if(f1.idx!=f2.idx){
    return (f1.idx<f2.idx);
  }else{
    return (f1.size>f2.size); 
  }
}

void CleanUpList(TString filelist, TString filename, bool print_uncut, TString FileName){

  ifstream input(filelist.Data());
  std::vector<filedata> file_vec; 
  filedata aux_data; 

  int name_sz = (int)filename.Sizeof();
  while(!input.eof()){
  
    TString line;
    input >> line;

    if(line.IsAlnum()){
      aux_data.size = line.Atoll(); 
      //std::cout << line.Data() << " ";
    }  

    if(line.Contains(filename.Data())){
      //std::cout << line.Data() << std::endl;
      TString str_aux = line;
      int idx_chop = line.Index("/store"); 
      line.Remove(0,idx_chop);   
      aux_data.name = line;
      int idx_aux = str_aux.Index(filename.Data());
      idx_aux+=name_sz;
      str_aux.Remove(0,idx_aux-1);
      str_aux.ReplaceAll(".root","");
      TString str_aux2 = str_aux;
      str_aux2.Remove(str_aux.Index("_"));
      str_aux.Remove(0,str_aux.Index("_")+1);

      aux_data.idx = str_aux2.Atoi();
      aux_data.subidx = str_aux.Atoi();
      if(print_uncut){
	//std::cout << aux_data.size << ", ";
	//std::cout << aux_data.name << ", ";
	//std::cout << aux_data.idx << ", ";
	//std::cout << aux_data.subidx << std::endl; 
      }

      file_vec.push_back(aux_data);  
    }
  }


  std::sort(file_vec.begin(), file_vec.end(), decision);

  int vec_sz = (int)file_vec.size();

  //if(print_uncut) std::cout << "*************** vec size after: " << (int)file_vec.size()  << " *************" << std::endl;  

  //ofstream myfile2;
  //myfile2.open ("prova.txt"); 
  if(print_uncut){
    for(int i=0; i<vec_sz;i++){
      //std::cout << file_vec.at(i).size << " " << file_vec.at(i).name.Data() <<std::endl;
      //myfile2   << file_vec.at(i).size << " " << file_vec.at(i).name.Data() <<std::endl;
    }
  }

  std::cout << std::endl; 
  std::cout << std::endl; 
  std::cout << "*************** Final List **************" << std::endl;  

  int no_repeat = 0;
  ofstream myfile;
  myfile.open (FileName); 
  for(int i=0; i<vec_sz;i++){
    if((i<vec_sz-1 && i>0 && file_vec.at(i).idx != file_vec.at(i-1).idx) || i==0 || i==vec_sz-1){
      //if((i<vec_sz-1 && i>0 && file_vec.at(i).idx != file_vec.at(i+1).idx && file_vec.at(i).idx != file_vec.at(i-1).idx) || (i==0 && file_vec.at(i).idx != file_vec.at(i+1).idx) || (i==vec_sz-1 && file_vec.at(i).idx != file_vec.at(i-1).idx)){
      //if(file_vec.at(i).size<1) std::cout << no_repeat << " " << file_vec.at(i).size << " " << file_vec.at(i).name.Data() <<std::endl;
      if(file_vec.at(i).size>1){
	myfile<<"'root://gridse2.pg.infn.it//storage/cms" << file_vec.at(i).name.Data() << "',"<<std::endl;
	no_repeat++;
      }
    }
  }
  
  int Nduplicates=0;
  std::cout << std::endl; 
  std::cout << std::endl; 
  std::cout << "*************** DOUBLE FILE ***************" << std::endl; 
  for(int i=0; i<vec_sz;i++){
    if(!((i<vec_sz-1 && i>0 && file_vec.at(i).idx != file_vec.at(i-1).idx) || i==0 || i==vec_sz-1)){
      //std::cout <<"'" << file_vec.at(i).name.Data() << "',"<<std::endl;
      Nduplicates++;
    }
  }
  if(Nduplicates==0) std::cout << "No duplicate file" << std::endl;


  std::cout << std::endl; 
  std::cout << std::endl; 
  std::cout << "*************** SUMMARY ***************" << std::endl; 
  std::cout << "Number of files without duplicates is " << no_repeat << std::endl;
  std::cout << "Total Number of files is " << (int)file_vec.size() << std::endl;
  //std::cout << "cnt: "<< cnt << std::endl;
  //std::cout << "Sum size: " << sum_size << std::endl;
  //std::cout << "Avg. size: " << (double)sum_size/cnt << std::endl;
}
