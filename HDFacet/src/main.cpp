#include "Decomposition.h"

int main(int argc, char* argv[])
{
	if(argc<2){
		cout<<"Argument list: "<<endl;
		cout<<"[1]: Path to model"<<endl;
		cout<<"[2]: Element size"<<endl;
		cout<<"[3]: Suffix"<<endl;
		cout<<"[4]: Density ratio"<<endl;
		cout<<"[5]: Step ratio"<<endl;
		cout<<"[6]: Steps"<<endl;
	}
	else{
		Decomposition decomtwo;

		decomtwo.initParameters(argv[3], stod(argv[4]), stod(argv[5]), stoi(argv[6]));
		decomtwo.buildMeshFromSurface(argv[1], stod(argv[2]));
		cout << "Triangulation Done ..." << endl;
		decomtwo.buildLaplacian();
		cout << "Laplacian Done ..." << endl;
		decomtwo.decompose();
		cout << "Decomposition Done ..." << endl;
		decomtwo.visualize();
		cout << "Visulization Done ..." << endl;
		decomtwo.integrate();
		decomtwo.write();
		cout << "Post Processing Done ..." << endl;

	}	

	return 1;
}
