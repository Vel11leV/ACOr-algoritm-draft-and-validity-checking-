//#include <thread>
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <cstdlib>
#include <ctime>
#include <cmath>


using namespace std;

#define pi 3,14159265358
/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	
	
	float Rosenbrock(int N, float *x)
	{
	double f =0;
	double X[N];
	for (int i=0;i<N;i++)
	X[i]=*(x+i);
	for (int i=0;i<(N-1);i++)
	f = f + 100*((X[i+1]-(X[i]*X[i]))*(X[i+1]-(X[i]*X[i])))+(X[i]-1)*(X[i]-1);
	
		return (float)f;
	}
	float Paraboloid (int N, float *x)
	{
	  float f =0;
	  for (int i=0;i<N;i++)
	  f = f + *(x+i)**(x+i);
		return f;
	}
	float Rastrigin(int N, float*x)
	{
		long double f = 0 ;
		long double X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		for(int i=0;i<N;i++)
		f = f + ((X[i]*X[i]+10)-10*cosl(((2*pi)*X[i])));
		
		return (float)f;
	}
	float Schaffer (int N, float*x)
	{
		float f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		f = 0.5+((sin((X[0]*X[0]+X[1]*X[1]))*sin((X[0]*X[0]+X[1]*X[1]))-0.5)/((1+0.001*(X[0]*X[0]+X[1]*X[1]))*(1+0.001*(X[0]*X[0]+X[1]*X[1]))));
		
		return f;
	}
	float Griewank (int N, float*x)
	{
		double f, fl=0, fr=1;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		for (int i=0;i<N;i++)
		{
		fl = fl + ((X[i])*X[i])/4000;
		fr = fr*(cos((X[i]/sqrt((i+1)))));
		}
	
		
		
		f = 1 + fl- fr;
		
		
		
		
		return  (float)(f);
	}
	float Himmelblau(int N, float*x)
	{
		double f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		f = (X[0]*X[0]+X[1]-11)*(X[0]*X[0]+X[1]-11)+(X[0]+ X[1]*X[1]-7)*(X[0]+ X[1]*X[1]-7);
		 
		return f;
	}
	
	float Martin_and_Gaddy(int N, float*x)
	{
		double f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		f = (X[0] - X[1])* (X[0] - X[1])+ ((X[0] + X[1]-10)/3)*((X[0] + X[1]-10)/3);
		
		return f;
	}
	
	
	/*
	float Easom(int N, float *x)
	{
		float f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		float p =-(X[0]-pi)*(X[0]-pi)-(X[1]-pi)*(X[1]-pi);
		f = (-cos(X[0]))*(cos(X[1]));//*exp(p);
		return f;
	}
	*/
	/*
	float Branin(int N, float *x)
	{
		long double f;
		long double X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		long double f1 = (X[1]-(5.1*X[0]*X[0]/(4*pi*pi))+(5/pi)*X[0]-6);
		f = f1*f1 + ((1-1/(8*pi))*cosl(X[0]))*10+10;
		
		
		return (float)f;
	}*/
	
	float Beale (int N, float *x)
	{
		double f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		float f1 = (1.5 -X[0]+X[0]*X[1]);
		float f2 = (2.25 -X[0]+X[0]*X[1]*X[1]);
		float f3 = (2.625 -X[0]+X[0]*X[1]*X[1]*X[1]);
		
		f = f1*f1 + f2*f2 + f3*f3;
		
		return (float)f;
	}
	
	float Three_hump_camel(int N, float *x)
	{
		float f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		f = 2*X[0]*X[0]-1.05*(X[0]*X[0])*(X[0]*X[0])+((X[0]*X[0])*((X[0]*X[0])*(X[0]*X[0]))/6)+X[0]*X[1]+X[1]*X[1];
		
		return f;
	}
	
	float Eggholder (int N, float *x)
	{
		float f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		float arg1 = sqrt(abs((X[0]/2) + (X[1]+47)));
		float arg2 = sqrt(abs((X[0] -(X[1]+47))));
		
		f = -(X[1]+47)*sin(arg1)-X[0]*sin(arg2);
		
		return (float)f;
	}
	
	float Bukin (int N, float *x)
	{
		float f;
		float X[N];
		for (int i=0;i<N;i++)
		X[i]=*(x+i);
		
		f = 100*sqrt(abs(X[1]-0.01*X[0]*X[0]))+ 0.01*abs(X[0]+10);
		
		return (float)f;
	}
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*-------------*/
class ACOR_algorithm
{
	public:
		
 	ACOR_algorithm(int Parameter_numb, int Archive_size, float Algor_param_q = 0.1 ,float Evaporation_ksi = 0.85, 
 				float optimum = 0.1,  unsigned int max_iterations = 1000,float (*Funct_to_opt) (int,float *) = NULL);
	~ACOR_algorithm ();
	void Param_bounds_setting(int param_number, float left, float right);
		int Solve ();
		
		
		/*-Random oscilators-*/
		float Rand_osc_uniform_0_1 ();
		float Rand_osc_uniform_a_b (float from_a, float to_b );
		int Rand_osc_roulette_wheel ();
		float Rand_osc_normal (float dispersion, float expected_value);
		/*-ACOrfunctions-*/
		float Standart_deviation (float evaporation, int sol_numb,int param_numb );
		float Solution_Weight (float q,int rank);
		void Probability_Solution_Weight ();

		void Sorting(void);
		bool New_solution_assembling ();
		void Initial_Archive_setting ();
		/*-Testfunction-*/
	
	
	private:
		float **theArchive ;//
		float **theParameter;//
		//float *theLimits;//&
		float *theNewSolution;
		
		float *theAlgor_param;//
		float *theOptimum_criteria;//
		int   *theParNumb;//
		int   *theArchiveSize;//
		unsigned int   *theMaxIterations;
		 int *theNumberFcalculation;
		bool  *theStatus;//
		float *theEvaporation;
		unsigned int *theTotalDoneIterationN;
		float (*theGoalFunction) (int, float *);
		float *theProbability_weight;

};
/*-------------*/
/*-------------*/
ACOR_algorithm::ACOR_algorithm(int Parameter_numb, int Archive_size,float Algor_param_q,float Evaporation_ksi, float optimum, unsigned int max_iterations,float (*FUNCT) (int,float*) )
				{
				theAlgor_param = new float;
				*theAlgor_param = Algor_param_q;
				
				theOptimum_criteria = new float;
				*theOptimum_criteria = optimum;
				
				theParNumb = new int;
				*theParNumb = Parameter_numb;
					
				theArchiveSize = new int;
				*theArchiveSize = Archive_size;
			
				theMaxIterations = new unsigned int;
				*theMaxIterations = max_iterations;
				
				theNumberFcalculation = new  int;
				*theNumberFcalculation  = 0;
				
				theStatus = new bool;
				*theStatus = false;
				
				theEvaporation = new float;
				*theEvaporation =Evaporation_ksi;
				
				theTotalDoneIterationN = new unsigned int;
				theTotalDoneIterationN = 0;
				
				theProbability_weight = new float[*theArchiveSize ];
				
				
				theNewSolution = new float [(*theParNumb)+2] ;
				for (int i=0;i<*theParNumb+2;i++)
				theNewSolution[i] = -1;
			
				theArchive = new float*[*theArchiveSize];
				for (int i = 0; i<*theArchiveSize ;i++)
				theArchive[i] = new float[*theParNumb+2];
			
				theParameter = new float* [*theParNumb];
				for (int i=0;i<*theParNumb; i++)
				theParameter[i] = new float[2];
				
				theGoalFunction =FUNCT;
				
				srand(time(NULL));
				}
/*----------*/
ACOR_algorithm::~ACOR_algorithm ()
		{
						
			for (int i=0;i<*theArchiveSize;i++)
					{
						delete [] theArchive [i];
						theArchive [i]= NULL;		
					}
			delete [] theArchive;
					theArchive = NULL ;
					
													
			for (int i=0;i<*theParNumb;i++)
				{
					delete [] theParameter [i];
					theParameter [i] = NULL;
				}
			delete [] theParameter;
					theParameter = NULL;
					
					
			delete 	[] theNewSolution;
					theNewSolution = NULL;
					
			delete [] theProbability_weight;	
					
			delete theAlgor_param;
				theAlgor_param = NULL;	
				
			delete theOptimum_criteria;
					theOptimum_criteria = NULL;
					
			delete theParNumb;
				theParNumb = NULL;		
							
			delete theArchiveSize;
				theArchiveSize =NULL;
				
			delete theMaxIterations;
				theMaxIterations = NULL;
				
			delete theNumberFcalculation;
				theNumberFcalculation = NULL;
						
			delete theStatus;
				theStatus = NULL;
				
			delete theEvaporation;
				theEvaporation =NULL;		
				
			delete theTotalDoneIterationN;
					theTotalDoneIterationN = NULL;
			
			cout<<"Destructor was called succesfully!\n";
			
		}
/*-----------*/
float ACOR_algorithm::Rand_osc_uniform_0_1 ()
{
	return	0.0001  *(rand() % 10001);
}
/*-----------*/
float ACOR_algorithm::Rand_osc_uniform_a_b (float from_a, float to_b )
{
	
	return from_a +(to_b-from_a)*ACOR_algorithm::Rand_osc_uniform_0_1 ();
}
/*-----------*/

float ACOR_algorithm::Rand_osc_normal (float dispersion, float expected_value)
{
	float S =-1;
	do{
	
	float rand1 = 0; 
	 rand1 = ACOR_algorithm::Rand_osc_uniform_a_b (-1, 1 );
	float rand2 = 0;
	 rand2 = ACOR_algorithm::Rand_osc_uniform_a_b (-1, 1 );
		
	S = rand1*rand1+rand2*rand2;
	if (S>0 && S<=1)
	{
	
	float sq = sqrt(-2*(log(S))/S);
	return dispersion*rand1*sq+expected_value;
	}
	else 
	S =-1;
	}while (S ==-1);

}
/*-----------*/
int ACOR_algorithm::Rand_osc_roulette_wheel ()
{
	
	int numb;
	float roulette[(*theArchiveSize)-1];
	float total_sum = 0;
	float roulette_arrow = 0;
	
		for (int i=0;i<*theArchiveSize;i++)
	       	total_sum = total_sum+theArchive[i][*theParNumb+1];
	    	       
		roulette [0]=theArchive[0][*theParNumb+1]/total_sum;
										
		for(int i=1;i<(*theArchiveSize)-1;i++)
				roulette [i]=(((theArchive[i][*theParNumb+1])/total_sum)+roulette [i-1]);
				 
		for(int i=0;i<*theArchiveSize-1;i++)
				roulette [i]=roulette [i]*100;
						
			roulette_arrow = Rand_osc_uniform_a_b (0, 100 );
		
		//for(int i=0;i<(*theArchiveSize-1);i++)
		//cout<<"\n"<<roulette [i] ;
		
	//	cout<<"\n"<<roulette_arrow ;
	if 	 (roulette_arrow <= roulette [0])
		return	numb = 0;
		
	if   (roulette_arrow>(roulette [(*theArchiveSize)-2]))
		return numb =  (*theArchiveSize)-1;
		
	for (int i=1;i<(*theArchiveSize)-1;i++)
		{
		if ((roulette_arrow <= roulette [i])&&(roulette_arrow > roulette [i-1]))
		return numb = i;
		}	
}

/*-----------*/
float ACOR_algorithm::Standart_deviation (float evaporation, int sol_numb, int param_numb)
{
	float dev = 0;
	float sum_dist = 0;
	for (int i=0;i<(*theArchiveSize);i++)
	{
		float dif = theArchive[i][param_numb]-theArchive[sol_numb][param_numb];
		
		if(dif <0)
		dif = -dif;
		
		sum_dist = sum_dist+dif;
	}
	
	dev = evaporation*sum_dist/(*theArchiveSize-1);
	if(dev <=0.00001)////////////////////////////////////////////////////////
	dev = 0.00001;///////////////////////////////////////////////////////////
/*	for (int i=0;i<(*theArchiveSize);i++)
	{
		if (i==sol_numb)
		continue;
		else
		{
		float distance = theArchive[i][param_numb]-theArchive[sol_numb][param_numb];
			if (distance<0)
				distance = -distance;
		dev = dev + distance/(*theArchiveSize-1);
		};
	}*/
	return dev;
}

/*----------*/
float ACOR_algorithm::Solution_Weight (float q,int rank)
{
	double weight = 0;
	double gfl = q;
	double rankfl = rank;
	double solnumb = *theArchiveSize;
	double part1 = 0, part2 = 0, part3 = 0;
	part1 = 1/(q*solnumb*sqrt(2*3.14159));
	part2 = -((rankfl-1)*(rankfl-1))/(2*q*q*solnumb*solnumb);
	part3 = exp(part2);
	weight = part3*part1;
	return (float)weight;
}
/*-----------*/
void ACOR_algorithm::Probability_Solution_Weight ()
{
	double sum = 0;
	for (int i=0;i<*theArchiveSize;i++)
	theProbability_weight[i]= Solution_Weight (*theAlgor_param,i+1);
	
	for (int i=0;i<*theArchiveSize;i++ )
	sum = sum + (double)theProbability_weight[i];
	
	for(int i=0;i<*theArchiveSize;i++)
	theProbability_weight[i]=(float)(theProbability_weight[i]/sum);
	
} 

/*-----------*/

	void ACOR_algorithm::Param_bounds_setting(int param_number, float left, float right)
	{
		theParameter[param_number-1][0] = left;
		theParameter[param_number-1][1] = right;
		
	}
	
/*-----------*/
void ACOR_algorithm::Sorting(void)
{
	 float temp [(*theParNumb)+1];
	for (int i = 0; i < (*theArchiveSize - 1); i++)
  {
    for (int j = (*theArchiveSize - 1); j > i; j--) // для всех элементов после i-ого
    {
      if (theArchive[j-1][*theParNumb] > theArchive[j][*theParNumb]) // если текущий элемент меньше предыдущего
      {
      	for (int k=0;k<(*theParNumb+1);k++)
        temp [k]=theArchive[j-1][k] ; // меняем их местами
        
        for (int k=0;k<(*theParNumb+1);k++)
        theArchive[j-1][k]=  theArchive[j][k]; 
        
        for (int k=0;k<(*theParNumb+1);k++)
        theArchive[j][k]=  temp [k]; 
      }
    }
  }
}
/*-----------*/

bool ACOR_algorithm::New_solution_assembling ()
{	
	for(int i=0;i<*theParNumb;i++)
		{			
				int Sol_Numb = Rand_osc_roulette_wheel ();
				
				float deviation = Standart_deviation (*theEvaporation, Sol_Numb,i);
				
				do	{		
					theNewSolution [i] = Rand_osc_normal ( deviation ,theArchive[Sol_Numb][i] );
					}	while (	theNewSolution [i]>theParameter[i][1]	||	theNewSolution [i]<theParameter[i][0]);											
		}	
		
		
	theNewSolution [*theParNumb] = (*theGoalFunction)(	*theParNumb, theNewSolution );
	(*theNumberFcalculation)++;
	return true;
}

/*-----------*/
void ACOR_algorithm::Initial_Archive_setting ()
{
	for (int i=0;i<*theArchiveSize;i++)
			{
		
			for (int j=0;j<(*theParNumb+1);j++ )
				{
						if(j<*theParNumb)
						{
						theNewSolution [j]=Rand_osc_uniform_a_b (theParameter[j][0], theParameter[j][1] );
						theArchive[i][j] =theNewSolution [j];
						}
						if(j==*theParNumb)
						{
						theNewSolution [j]=(*theGoalFunction)(	*theParNumb, theNewSolution );
						theArchive[i][j] = theNewSolution [j];
						}
				}
			}
		Sorting();	
		for (int k=0;k<*theArchiveSize;k++)
		theArchive[k][*theParNumb+1] = Solution_Weight (*theAlgor_param,k+1);
		Probability_Solution_Weight ();
		/*
		for (int j=0;j<*theArchiveSize;j++)
		{
		
		for (int i=0;i<*theParNumb+2;i++)
		cout<<setprecision(3)<<theArchive[j][i]<<"\t"<<" ";
		cout<<"\n";
		}
		cout<<"\n\n\n";
		for(int i=0;i<*theParNumb+2;i++)
		cout<<theNewSolution [i]<<"  ";*/
}
/*-----------*/



int ACOR_algorithm::Solve ()

{
	
		Initial_Archive_setting ();
		*theNumberFcalculation = 0;
		Probability_Solution_Weight ();
		float arch[2][*theParNumb+2];
		do
		{
			
	
			
		do{
			New_solution_assembling ();
		
			if (*theNumberFcalculation>=*theMaxIterations)
			break;
		} while ((theArchive[(*theArchiveSize)-2][*theParNumb])<(theNewSolution [*theParNumb]));
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[0][i] = theNewSolution[i];
		
		
		do{
			New_solution_assembling ();
			if (*theNumberFcalculation>=*theMaxIterations)
			break;
		} while ((theArchive[(*theArchiveSize)-2][*theParNumb])<(theNewSolution [*theParNumb]));
		
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[1][i] = theNewSolution[i];
		
		
		/*
		do{
			New_solution_assembling ();
		
		} while ((theArchive[(*theArchiveSize)-5][*theParNumb])<(theNewSolution [*theParNumb]));
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[2][i] = theNewSolution[i];
		do{
			New_solution_assembling ();
		
		} while ((theArchive[(*theArchiveSize)-5][*theParNumb])<(theNewSolution [*theParNumb]));
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[3][i] = theNewSolution[i];
		do{
			New_solution_assembling ();
		
		} while ((theArchive[(*theArchiveSize)-5][*theParNumb])<(theNewSolution [*theParNumb]));
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[4][i] = theNewSolution[i];
		
		do{
			New_solution_assembling ();
		
		} while ((theArchive[(*theArchiveSize)-5][*theParNumb])<(theNewSolution [*theParNumb]));
		
		for (int i=0;i<(*theParNumb)+1;i++)
		arch[4][i] = theNewSolution[i];
		*/
		
		
		for (int i=0;i<(*theParNumb)+1;i++)
		{
		
		theArchive[(*theArchiveSize)-1][i] =arch[0][i];
		theArchive[(*theArchiveSize)-2][i] =arch[1][i];
		/*theArchive[(*theArchiveSize)-3][i] =arch[2][i];
		theArchive[(*theArchiveSize)-4][i] =arch[3][i];
		theArchive[(*theArchiveSize)-5][i] =arch[4][i];*/
		}
		
						Sorting();
						cout<<"\n"<<theArchive[0][*theParNumb]<<"\n";
						for(int i=0;i<(*theParNumb+1);i++)
							{
							cout<<theArchive[0][i]<<"   ";
							}
				if (*theNumberFcalculation>=*theMaxIterations)
				{
				*theNumberFcalculation = (-1);
				break;
				}
		} while (theArchive[0][*theParNumb]>*theOptimum_criteria );
		
		cout <<"\nOptimal parametres values are : \n";
		for(int i=0;i<(*theParNumb);i++)
							{
							cout<<theArchive[0][i]<<"   ";
							}
		cout<<"\n The best GF value is : "<<theArchive[0][*theParNumb]<<"\n";
	/*
		for (int i=0;i<*theArchiveSize;i++)
		{
			for (int j=0;j<(*theParNumb)+2;j++)
			cout<<setprecision(1)<<theArchive[i][j]<<"\t";
			cout<<setprecision(2)<<theProbability_weight[i]<<"\n";
		
		}*/
		/*
		fstream fout("cppstudio.txt"); 
			
		fout<<"_________________________________________"<<'\n';					
		fout<<"OPTIMAL PARAMETERS ARE:"<<'\n';	
		fout<<"_________________________________________"<<'\n';	
		for(int i=0;i<(*theParNumb);i++)
		fout<<theArchive[0][i]<<"  \n";
		fout<<"_________________________________________"<<'\n';	
		fout<<"THE GOAL FUNCTION VALUE IS:"<<'\n';	
		fout<<"_________________________________________"<<'\n';	
		fout<<theArchive[0][*theParNumb]<<"  \n";
		fout<<"_________________________________________"<<'\n';	
		
		
		/*cout<<'\n';	
		for (int i=0;i<(*theParNumb)+1;i++)
			{
			
			cout<<theNewSolution [i]<<"  ";
			}
			/*
			float d = 0, g = 0;
		
		cout<<"\n";	
		for (int i=0;i< 5; i++)
		{
			d = d + theProbability_weight[i];
		}
		cout<<"\n"<<"First 5 :"<<d<<"\n";
		
		for (int i=0;i< 10; i++)
		{
			g = g + theProbability_weight[i];
		}
		cout<<"First 10 :"<<g<<"\n\n";
		
		*/
		/*
		fout<<"___________________________________________"<<'\n';	
		fout<<"THE NUMBER OF GOAL FUNCTION CALCULATIONS IS:"<<'\n';	
		fout<<"___________________________________________"<<'\n';	
		fout<<*theNumberFcalculation<<"\n";
		fout<<"___________________________________________"<<'\n';	
		
		fout.close(); // закрываем файл
		*/		
	return *theNumberFcalculation;
}




int main()
{
	float q = 0.1, ksi = 0.85;
	float archive_size = 50;
	float threshhold=100;
	float bound_l = -30, bound_r = 30;
	int par_numb = 30;
	int Iterations = 5;
	float x[2]={100,100};
	char par = NULL;
	float (*Funct) (int,float *) = Rosenbrock;
	int F_opt = 0;
	string  Function = "Rosenbrock";
	char anyk;

	//Alg.Solve();
	/*
	ofstream fout("cppstudio.txt",ios_base::app); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
    fout << "\nПривет мир!"; // запись строки в файл
    fout.close(); // закрываем файл
    fout.open("cppstudio.txt", ios_base::app|ios_base::ate);
    fout << "\nЗапись в файл "; 
    fout.close();*/

	
	do	{
		char in_ch = 'a';
		cout<<"===================================================================\n";
		cout<<"===================================================================";
		cout<<"\n\nStandart algorithm parametres are:\nq = "<<q<< "; ksi = "<< ksi<< "; archive size = "<< archive_size<<". Two ants per iteration. \n";
		cout<<"\nInitial parametres setting press: [c] + Enter; to ignor - [any key] + Enter. \n";
		cin>>in_ch;
		if (in_ch=='c')
		{
			cout<<"Enter the q [0.001,1]: ";
			cin>>q;
			cout<<"Enter the ksi [0.001,100]: ";
			cin>>ksi;
			cout<<"Enter the archive size [2...100]: ";
			cin>>archive_size;
			cout<<"New parametres were suxcessfully set...\n";
		}
		else
			cout<< "Standart parametres are in use...\n";
		
		
		cout<<"\n[0]-Rosenbrock; [1]-Paraboloid; [2]-Rastrigin; [3]-Schaffer; [4]-Griewank; [5]-Himmelblau; [6]-Martin_and_Gaddy; [7]-Beale ; [8]-Eggholder; [9]-Bukin	\n";
		cout<<"Enter the function to be optimized: ";
		cin>>F_opt;
		
		cout<<"Enter the number of alg. iterations: ";
		cin>>Iterations;
		//cout<<"\n";
		
		cout<<"Enter the GF threshhold: ";
		cin>>threshhold;
		//cout<<"\n";
		
		cout<<"Enter the number of parametres: ";
		cin>>par_numb;
		//cout<<"\n";
		
		cout<<"Enter the left side parametre limitation ( for instance : -1.1):  ";
		cin>>bound_l;
		//cout<<"\n";
		
		cout<<"Enter the right side parametre limitation (for instance : 1.1):  ";
		cin>>bound_r;
		//cout<<"\n";
		cout<<"Parametres were set correctly. To start calculation press any key...";
		cin>>anyk;
		cout<<"\n";
		
		
			switch (F_opt) 
				{
					case 0: Funct = Rosenbrock;
							Function = "Rosenbrock";
							break;
					case 1: Funct = Paraboloid;
							Function = "Paraboloid";
							break;
					case 2: Funct = Rastrigin;
							Function = "Rastrigin";
							break;	
					case 3: Funct = Schaffer;
							Function = "Schaffer";
							break;
					case 4: Funct = Griewank;
							Function = "Griewank";
							break;	
					case 5: Funct = Himmelblau;
							Function = "Himmelblau";
							break;
					case 6: Funct = Martin_and_Gaddy;
							Function = "Martin_and_Gaddy";
							break;
					case 7: Funct = Beale;
							Function = "Beale";
							break;	
					case 8: Funct = Eggholder;
							Function = "Eggholder";
							break;
					case 9: Funct = Bukin;
							Function = "Bukin";
							break;		
				}
		 int N[Iterations];
		ACOR_algorithm Alg(par_numb ,archive_size, q ,ksi, threshhold,500000, Funct);
	
		for (int i=1;i<par_numb+1;i++)
		Alg.Param_bounds_setting(i, bound_l , bound_r);
		
		Alg.Initial_Archive_setting ();
		
		for (int i=0;i<Iterations;i++)
		N[i] = Alg.Solve();
		cout<<"\n____________________________________________\n";
		cout<<"Algorithm parametres: q = "<<q<<"; ksi = "<<ksi<<"; archive size = "<<archive_size<<". Two ants per iteration.\n";
		cout<<"\nFUNCTION:  "<<Function<<"\n";
		cout<<"\nNumber of parametres is: "<< par_numb<<"\n";
		cout<<"\nParametre value is: ["<<bound_l<<" , "<<bound_r <<"]\n";
		cout<<"\nThreshhold is: "<<threshhold<<"\n";
		for (int i=0;i<Iterations;i++)
		cout<<"\nIteration "<<(i+1)<<" = "<<N [i];
		
		int SUM =0;
		int dev =0;
		for (int i=0;i<Iterations;i++)
			{
			if (N[i]>0)
			{
			 SUM = SUM + N[i];
			  dev ++;
			}
			else continue;
			}
		if (dev==0)
		dev = 1;
	    cout<<"\nThe average number of GF calculations is: "<< SUM/dev<<"\n";
	    
	    								char save = 'a';
	    								cout<<"\nTo save press [s+Enter], to ignor press[any key + Enter]: \n";
	    								cin>> save;
	    								if (save=='s')
	    								{
	    								
	    								fstream fout("ACOr_Result.txt",ios_base::out|ios_base::app); 
	    								fout<<"\n____________________________________________\n";
	    								fout<<"\nAlgorithm parametres: q = "<<q<<"; ksi = "<<ksi<<"; archive size = "<<archive_size<<". Two ants per iteration.\n";
										fout<<"\nFUNCTION:  "<<Function<<"\n";
										fout<<"\nNumber of parametres is: "<< par_numb<<"\n";
										fout<<"\nParametre value is: ["<<bound_l<<" , "<<bound_r <<"]\n";
										fout<<"\nThreshhold is: "<<threshhold<<"\n";
										for (int i=0;i<Iterations;i++)
										fout<<"\nIteration "<<(i+1)<<" = "<<N [i];
																		
									    fout<<"\nThe average number of GF calculations is: "<< SUM/dev<<"\n";
									    
	    								fout.close();
	    								cout<< "\n The data was saved sucsesfully....\n";
										}
										else
										cout<<"Saving ignored...";
		cout<<"\nTo quit the program press -[q] + Enter, to continue - press  any key + Enter...";
		cin>>par;
		
		
		
	}while (par != 'q'); 
	
	
			
	return 0;
}





