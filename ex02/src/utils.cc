#include "utils.h"
double delta;
using namespace std;

//////////////////////////////////////////////////
void MatrixPrint(vector< map<int, double> > &A, vector< map<int, double> > &M, const int NX){

     map<int,double>::iterator it;
     std::ofstream myfile;
     myfile.open("../output/A.txt");

     for(int i=0;i<NX;i++){
        for (it=A[i].begin(); it!=A[i].end(); ++it){
            myfile<< i <<" "<< it->first << " " << it->second< <"\n";
        }
     }

     myfile.close();

     myfile.open("../output/M.txt");

     for(int i=0;i<NX;i++){
        for (it=M[i].begin(); it!=M[i].end(); ++it){
            myfile<< i <<" "<< it->first << " " << it->second <<"\n";
        }
     }
     myfile.close();
}

//////////////////////////////////////////////////
void initialize(double* vect, const int memSize, const double value){ 

    for(int i=0; i<memSize; i++){   
        vect[i] = value;
    }

}

//////////////////////////////////////////////////
//Computes L2 norm
double L2N(double * const resid, const int memSize){ 

    double sum = 0.0;

    for(int i=0; i<memSize; i++){
        sum += resid[i]*resid[i];
    }

    return sqrt(sum);
}


////////////////////////////////////////////////////
bool ReadFile(const std::string &name, double* &vertex, int & nv, int* &face, int & nf){

    std::ifstream input( name.c_str() );
    std::string line;
    int state=0;
    int nr;
    
    if (input!=NULL) state=1;//Check if the file is good
    
    for (int i=0;i<2;i++){
    	getline( input, line);
    	istringstream issnr(line);
    	std::string sub;
        issnr >>std::ws>>sub>>std::ws;

    	if(0==i){
            nv=atoi(sub.c_str());
    		nr=nv;
            allocMem(&vertex , nv*2);

            // Init vertexes to zero
            for(int t=0; t<nv*2; t++) 
                vertex[t]=0.0;
	    }
	    else{nf=atoi(sub.c_str());
	        nr=nf;
             
        allocMem(& face, nf*3);

        // Init faces to zero
        for(int p=0; p<nf*3; p++) 
            face[p]=0.0;
	}

	getline( input, line);
	
    for(int j=0;j<nr;j++){	
	   
        getline( input, line);
        istringstream iss(line);
        std::string sub1, sub2, sub3;
        iss >> std::ws >> sub1 >> std::ws
                >> sub2 >> std::ws >> std::ws >>sub3;

        if(0==i){
            double x=atof(sub2.c_str());
            double y=atof(sub3.c_str());
            vertex[j*2]=x;
            vertex[j*2+1]=y;
        }
	    else{
            int index0=atoi(sub1.c_str());
    		int index1=atoi(sub2.c_str());
    		int index2=atoi(sub3.c_str());
    		face[j*3]=index0;
    		face[j*3+1]=index1;
    		face[j*3+2]=index2;}				
	    }
    }
     
    if (state==1 )return true;
    else return false;
		
}

////////////////////////////////////////////////////
void faceprint(int * const array, const int n){

    for(int i=0;i<n;i++){
        for(int j=0;j<3;j++){
            cout<< array[i*3+j] << " ";
        }
        cout<<endl;
     }
}

////////////////////////////////////////////////////
void vertexprint(double * const array, const int n){

     int k=0;
     for(int i=0;i<n;i++){
        cout< <k++ << " ";
        for(int j=0;j<2;j++){
            cout<< array[i*2+j ]<< " ";
        }
        cout<<endl;
     }
}

///////////////////////////////////////////////////
void solvek(double * k, double * const array, const int n,const double delt){

     for(int i=0;i<n;i++){
        k[i]=(100+delt)*pow(M_E,-50*(array[i*2]*array[i*2]+array[i*2+1]*array[i*2+1]))-100;
     }
     delta=delt;
}

//////////////////////////////////////////////////
void kprint(double * const k ,double * const array, const int n){

     std::ofstream myfile;
     myfile.open("../output/ksq.txt");
     for(int j=0;j<n;j++){
         myfile<<array[j*2]<<" "<<array[j*2+1]<<" "<<k[j]<<"\n";
     }
     
}

//////////////////////////////////////////////////
void CreateNeighbourList(int * nneighb, int * &neighb, int * const f, const int nv, const int nf, int & nn){

    //we first count the number of neighbours in order to allocate memory for the array
   for(int i=0;i<nv;i++){ 

        int n=0;
        for(int j=0;j<nf;j++){
            if(i==f[j*3]||i==f[j*3+1]||i==f[j*3+2]) n++;
        }

	    nneighb[2*i]=n;  //the number of neighbours for that vertex
        nn+=n;
    }
     
    allocMem(& neighb,nn);  //the size of neighb is the sum of the neighbours of all vertexes
    int nr=0;
     
    for(int i=0;i<nv;i++){
	
        nneighb[2*i+1]=nr;  //the position of the first neighbour of the vertex in the neighb array 
        
        for(int j=0;j<nf;j++){
            if(i==f[j*3]||i==f[j*3+1]||i==f[j*3+2]) neighb[nr++]=j;
        }
    }  
}

///////////////////////////////////////////////////
void getCoord(double * const vertex, const int poz, double& x, double& y){
     
     x=vertex[poz*2];
     y=vertex[poz*2+1];

}

///////////////////////////////////////////////////
void getVertexes(int * const face, const int nr, int& v1, int& v2, int& v3){
    
    v1=face[nr*3];
    v2=face[nr*3+1];
    v3=face[nr*3+2];

}

//////////////////////////////////////////////////
void getCorners(double * corners,double * vertex,int * face, const int pos){
    
    int v1=0;
    int v2=0;
    int v3=0;
    double x=0;
    double y=0;
    getVertexes(face, pos, v1, v2, v3);

    getCoord(vertex, v1, x, y);
    corners[0]=x;
    corners[1]=y;
    getCoord(vertex, v2, x, y);
    corners[2]=x;
    corners[3]=y;
    getCoord(vertex, v3, x, y);
    corners[4]=x;
    corners[5]=y;

}

///////////////////////////////////////////////////
double comGrad(const double x, const double y){
    
    return  (   (100+delta) *  pow( M_E,-50*(x*x+y*y) ) - 100   );

}

//////////////////////////////////////////////////
void copyArray(double * const source, double * const dest, const int n){
   
    for(int i=0;i<n;i++)
        dest[i]=source[i];

}

//////////////////////////////////////////////////
void eigprint(double * const eig ,double * const array, const int n){

     std::ofstream myfile;
     myfile.open("../output/eigenmode.txt");
     
     for(int j=0;j<n;j++){
         myfile<<array[j*2]<<" "<<array[j*2+1]<<" "<<eig[j]<<"\n";
     }
     
}
