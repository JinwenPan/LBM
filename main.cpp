#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <any>
using namespace std;
const double Wc = (double)4 / 9;
const double Wn = (double)1 / 9; const double Ws = (double)1 / 9; const double Ww = (double)1 / 9; const double We = (double)1 / 9;
const double Wnw = (double)1 / 36; const double Wne = (double)1 / 36; const double Wsw = (double)1 / 36; const double Wse = (double)1 / 36;

namespace lbm {
	using uint = unsigned int;

	template <typename Type, uint Cellsize>
	class Grid {
	public:
		inline Grid() = default;
		inline Grid(uint xsize, uint ysize);
		inline ~Grid() = default;
		inline Type& operator()(uint x, uint y, uint f);
		inline Type operator()(uint x, uint y, uint f) const;
		void swap(Grid& grid);
	private:
		uint xsize_;
		uint ysize_;
		vector<Type> data_;
	};

	template <typename Type, uint Cellsize>
	Grid<Type, Cellsize>::Grid(uint xsize, uint ysize)
		: xsize_(xsize)
		, ysize_(ysize)
		, data_(vector<Type>(Cellsize*xsize*ysize))
	{}

	template <typename Type, uint Cellsize>
	inline Type& Grid<Type, Cellsize>::operator()(uint x, uint y, uint f) {
		assert(x < xsize_ && y < ysize_ && f < Cellsize);
		return data_[y * xsize_ * Cellsize + x * Cellsize + f];
	}

	template <typename Type, uint Cellsize>
	inline Type Grid<Type, Cellsize>::operator()(uint x, uint y, uint f) const{
		assert(x < xsize_ && y < ysize_&& f < Cellsize);
		return data_[y * xsize_ * Cellsize + x * Cellsize + f];
	}

	using PDF_Field = Grid<double, 9>;
	using Velocity_Field = Grid<double, 2>;
	using Density_Field = Grid<double, 1>;
	using Flag_Field = Grid<uint, 1>;

	template <typename Type, uint Cellsize>
	void Grid<Type, Cellsize>::swap(Grid& grid) {
		using std::swap;
		swap(xsize_, grid.xsize_);
		swap(ysize_, grid.ysize_);
		swap(data_, grid.data_);
	}


	template < typename Type , uint Cellsize >
	inline void swap ( Grid <Type , Cellsize >& a,Grid <Type , Cellsize >& b ) {
	a. swap ( b );
	}

//===============================================Flag initialization================================================
	void Flag_Initialization(Flag_Field& grid_flag, uint sizex, uint sizey, uint spherex, uint spherey, uint diameter) {
		for (uint i = 0; i < sizex; ++i) {
			grid_flag(i, 0, 0) = 1;                    //setting flag 1 for non slip upper and lower boundaries
			grid_flag(i, sizey - 1, 0) = 1;
		}
		for (uint j = 1; j < sizey - 1; ++j) {
			grid_flag(0, j, 0) = 2;                     //setting flag 2 & 3 for left and the righ wall
			grid_flag(sizex - 1, j, 0) = 3;
		}
		for (uint j = 1; j < sizey - 1; ++j) {
			for (uint i = 1; i < sizex - 1; ++i) {      //setting flag zero for fluid cells
				grid_flag(i, j, 0) = 0;
			}
		}
		float radius = (float)diameter / 2;             //setting flag 1 for the obstacle
		if ((diameter % 2) != 1) {
			for (uint j = 0; j < diameter; ++j) {
				for (uint i = 0; i < diameter; ++i) {
					if ( pow(0.5 + i - radius, 2) + pow(0.5 + j - radius, 2) < pow(radius, 2) ) {
						grid_flag(spherex - radius + 1 + i, spherey - radius + 1 + j, 0) = 1;
					}
				}
			}
		}
		else {
			for (uint j = 0; j < (diameter + 1); ++j) {
				for (uint i = 0; i < (diameter + 1); ++i) {
					if ( pow(i - radius, 2) + pow(j - radius, 2) < pow(radius, 2) ) {
						grid_flag(spherex - radius + 0.5 + i, spherey - radius + 0.5 + j, 0) = 1;
					}
				}
			}
		}
	}
//===============================================PDF initialization================================================
	void PDF_Initialization(PDF_Field& grid_pdf, uint sizex, uint sizey, uint spherex, uint spherey, uint diameter) {
		for (uint i = 0; i < sizex; ++i) {
			for (uint k = 0; k < 9; ++k) {                        // Initializing the PDFs for the helper cells at the upper and lower boundary  (setting to zero)
				grid_pdf(i, 0, k) = 0;
				grid_pdf(i, sizey - 1, k) = 0;
			}
		}
		for (uint j = 1; j < sizey - 1; ++j) {
			for (uint k = 0; k < 9; ++k) {                       //Initializing the PDFs for the helper cells at the left and right bounadry
				grid_pdf(0, j, k) = 0;
				grid_pdf(sizex - 1, j, k) = 0;
			}
		}
		for (uint j = 1; j < sizey - 1; ++j) {
			for (uint i = 1; i < sizex - 1; ++i) {                //Initializing PDFs for fluid cells
				grid_pdf(i, j, 0) = Wc;
				grid_pdf(i, j, 1) = Wn; grid_pdf(i, j, 2) = Ws; grid_pdf(i, j, 3) = Ww; grid_pdf(i, j, 4) = We;
				grid_pdf(i, j, 5) = Wnw; grid_pdf(i, j, 6) = Wne; grid_pdf(i, j, 7) = Wsw; grid_pdf(i, j, 8) = Wse;
			}
		}
		float radius = (float)diameter / 2;                      //Initializing PDFs for the obstacle (setting to zero)
		if ((diameter % 2) != 1) {
			for (uint j = 0; j < diameter; ++j) {
				for (uint i = 0; i < diameter; ++i) {
					if ( pow(0.5 + i - radius, 2) + pow(0.5 + j - radius, 2) < pow(radius, 2) ) {
						for (uint k = 0; k < 9; ++k) {
							grid_pdf(spherex - radius + 1 + i, spherey - radius + 1 + j, k) = 0;
						}
					}
				}
			}
		}
		else {
			for (uint j = 0; j < (diameter + 1); ++j) {
				for (uint i = 0; i < (diameter + 1); ++i) {
					if (pow(i - radius, 2) + pow(j - radius, 2) < pow(radius, 2)) {
						for (uint k = 0; k < 9; ++k) {
							grid_pdf(spherex - radius + 0.5 + i, spherey - radius + 0.5 + j, k) = 0;
						}
					}
				}
			}
		}
	}
//===============================================Density================================================
	void Density(Density_Field& grid_density, PDF_Field& grid_pdf,Flag_Field& grid_flag, uint sizex, uint sizey) {
		for (uint j = 0; j < sizey; ++j) {
			for (uint i = 0; i < sizex; ++i) {              

				double tmp = 0;
				for (uint k = 0; k < 9; ++k) {                    //Calculating density for fluid cells
					if (grid_flag(i, j, 0) == 0) {
					tmp += grid_pdf(i, j, k);
					}
				grid_density(i, j, 0) = tmp;                            
				}

			}
		}
	}
//===============================================Velocity===================================================
	void Velocity(Velocity_Field& grid_velocity, PDF_Field& grid_pdf,Flag_Field& grid_flag, uint sizex, uint sizey) {
		for (uint j = 0; j < sizey; ++j) {
			for (uint i = 0; i < sizex; ++i) {
				if (grid_flag(i, j, 0) == 0) {
					double tmpx = 0;
					tmpx += grid_pdf(i, j, 4);
					tmpx += grid_pdf(i, j, 6);                         //Calculating velocity for fluid cells
					tmpx += grid_pdf(i, j, 8);
					tmpx -= grid_pdf(i, j, 3);
					tmpx -= grid_pdf(i, j, 5);
					tmpx -= grid_pdf(i, j, 7);
					grid_velocity(i, j, 0) = tmpx;
					double tmpy = 0;
					tmpy += grid_pdf(i, j, 1);
					tmpy += grid_pdf(i, j, 5);
					tmpy += grid_pdf(i, j, 6);
					tmpy -= grid_pdf(i, j, 2);
					tmpy -= grid_pdf(i, j, 7);
					tmpy -= grid_pdf(i, j, 8);
					grid_velocity(i, j, 1) = tmpy;
				}

			}
		}
	}
//===============================================PDF Equilibrium================================================
	void PDF_Equilibrium(PDF_Field& grid_pdf_eq, Density_Field& grid_density, Velocity_Field& grid_velocity, Flag_Field& grid_flag, uint sizex, uint sizey) {
		for (uint j = 1; j < sizey-1; ++j) {
			for (uint i = 1; i < sizex-1; ++i) {
				if (grid_flag(i, j, 0) == 0) {        //Calculating PDF equilibrium for fluid cells 
					double p = grid_density(i, j, 0);
					double ux = grid_velocity(i, j, 0);
					double uy = grid_velocity(i, j, 1);
					double usquared = grid_velocity(i, j, 0) * grid_velocity(i, j, 0) + grid_velocity(i, j, 1) * grid_velocity(i, j, 1);
					grid_pdf_eq(i, j, 0) = Wc * (p - 1.5 * usquared);
					grid_pdf_eq(i, j, 1) = Wn * (p + 3 * uy + 4.5 * uy * uy - 1.5 * usquared);
					grid_pdf_eq(i, j, 2) = Ws * (p - 3 * uy + 4.5 * uy * uy - 1.5 * usquared);
					grid_pdf_eq(i, j, 3) = Ww * (p - 3 * ux + 4.5 * ux * ux - 1.5 * usquared);
					grid_pdf_eq(i, j, 4) = We * (p + 3 * ux + 4.5 * ux * ux - 1.5 * usquared);
					grid_pdf_eq(i, j, 5) = Wnw * (p + 3 * (-ux + uy) + 4.5 * pow(-ux + uy, 2) - 1.5 * usquared);
					grid_pdf_eq(i, j, 6) = Wne * (p + 3 * (ux + uy) + 4.5 * pow(ux + uy, 2) - 1.5 * usquared);
					grid_pdf_eq(i, j, 7) = Wsw * (p + 3 * (-ux - uy) + 4.5 * pow(-ux - uy, 2) - 1.5 * usquared);
					grid_pdf_eq(i, j, 8) = Wse * (p + 3 * (ux - uy) + 4.5 * pow(ux - uy, 2) - 1.5 * usquared);
				}
			}
		}
	}	
//===============================================Collide===========================================================
	void Collide(PDF_Field& grid_pdf, PDF_Field& grid_pdf_eq, Flag_Field& grid_flag, double t, uint sizex, uint sizey) {
		for (uint j = 1; j < sizey - 1; ++j) {
			for (uint i = 1; i < sizex - 1; ++i) {                         //Performing collide step for fluid cells
				if (grid_flag(i, j, 0) == 0) {
					for (uint k = 0; k < 9; ++k) {
						grid_pdf(i, j, k) = grid_pdf(i, j, k) - grid_pdf(i, j, k) / t + grid_pdf_eq(i, j, k) / t;
					}
				}
			}
		}
	}
//===============================================Boundary Handling============================================================
	void Boundary_handling(PDF_Field& grid_pdf, Flag_Field& grid_flag, Velocity_Field& grid_velocity, double uin, uint sizex, uint sizey, uint spherex, uint spherey, uint diameter){

	//===============================================North Wall===========================================================
		grid_pdf(0,sizey-1,8)=grid_pdf(1,sizey-2,5);              
		grid_pdf(sizex-1,sizey-1,7)=grid_pdf(sizex-2,sizey-2,6); 

		for (uint i=1; i<sizex-1; ++i){                                     //Sending fluid cells PDF's to the helper cells at the upper boundary
			grid_pdf(i,sizey-1,2)=grid_pdf(i,sizey-2,1);
			grid_pdf(i,sizey-1,7)=grid_pdf(i-1,sizey-2,6);
			grid_pdf(i,sizey-1,8)=grid_pdf(i+1,sizey-2,5);

		}

	//===============================================South Wall===========================================================
		grid_pdf(0,0,6)=grid_pdf(1,1,7);    
		grid_pdf(sizex-1,0,5)=grid_pdf(sizex-2,1,8);                //Sending fluid cells PDF's to the helper cells at the lower boundary

		for (uint i=1; i<sizex-1; ++i){    
			grid_pdf(i,0,1)=grid_pdf(i,1,2);
			grid_pdf(i,0,6)=grid_pdf(i+1,1,7);
			grid_pdf(i,0,5)=grid_pdf(i-1,1,8);
		}
	//===============================================West Wall===========================================================
		for (uint j=1; j<sizey-1; ++j){
			grid_pdf(0,j,4)=grid_pdf(1,j,3)+6*Ww*uin;                //Sending fluid cells PDF's to the helper cells (Velocity boundary)
			grid_pdf(0,j,6)=grid_pdf(1,j+1,7)+6*Wsw*uin;
			grid_pdf(0,j,8)=grid_pdf(1,j-1,5)+6*Wne*uin;
		}


	//===============================================East Wall===========================================================
		for (uint j=1; j<sizey-1;  ++j){
		grid_pdf(sizex-1,j,3)=-grid_pdf(sizex-2,j,4)+2*We*(1+4.5*(pow(grid_velocity(sizex-2,j,0),2)) - 1.5* ((pow(grid_velocity(sizex-2,j,0),2))+(pow(grid_velocity(sizex-2,j,1),2))));
		grid_pdf(sizex-1,j,5)=-grid_pdf(sizex-2,j+1,8)+2*Wse*(1+4.5*(pow(grid_velocity(sizex-2,j+1,0)-grid_velocity(sizex-2,j+1,1),2))- 1.5* ((pow(grid_velocity(sizex-2,j+1,0),2))+(pow(grid_velocity(sizex-2,j+1,1),2))));
		grid_pdf(sizex-1,j,7)=-grid_pdf(sizex-2,j-1,6)+2*Wne*(1+4.5*(pow(grid_velocity(sizex-2,j-1,0)+grid_velocity(sizex-2,j-1,1),2))- 1.5* ((pow(grid_velocity(sizex-2,j-1,0),2))+(pow(grid_velocity(sizex-2,j-1,1),2))));
		}                                                             //Sending fluid cells PDF's to the helper cells (Denisty boundary)

	//===============================================Obstacle===========================================================
		uint sphere_bottom = spherey - diameter/2;
		uint sphere_top = spherey + diameter/2;
		uint sphere_left = spherex - diameter/2;
		uint sphere_right = spherex + diameter/2;
		for (uint j=sphere_bottom-2; j<sphere_top+2; ++j){          //Sending fuild cells PDF's to the adjacent obstacle cells
			for (uint i=sphere_left-2; i<sphere_right+2; ++i){
				if (grid_flag(i,j,0)==1){
					if(grid_flag(i+1,j,0)==0){
						grid_pdf(i,j,4)=grid_pdf(i+1,j,3);
					}
					if(grid_flag(i+1,j+1,0)==0){
						grid_pdf(i,j,6)=grid_pdf(i+1,j+1,7);
					}
					if(grid_flag(i+1,j-1,0)==0){
						grid_pdf(i,j,8)=grid_pdf(i+1,j-1,5);
					}
					if(grid_flag(i-1,j,0)==0){
						grid_pdf(i,j,3)=grid_pdf(i-1,j,4);
					}
					if(grid_flag(i-1,j+1,0)==0){
						grid_pdf(i,j,5)=grid_pdf(i-1,j+1,8);
					}
					if(grid_flag(i-1,j-1,0)==0){
						grid_pdf(i,j,7)=grid_pdf(i-1,j-1,6);
					}
					if(grid_flag(i,j+1,0)==0){
						grid_pdf(i,j,1)=grid_pdf(i,j+1,2);
					}
					if(grid_flag(i,j-1,0)==0){
						grid_pdf(i,j,2)=grid_pdf(i,j-1,1);
					}
				}
			}
		}


	}

//===============================================Streaming Pull============================================================
	void Pulling(PDF_Field& grid_pdf_dest,PDF_Field& grid_pdf, Flag_Field& grid_flag, uint sizex, uint sizey){
		for(uint j=1; j < sizey-1; ++j){
			for (uint i = 1; i <sizex-1 ; ++i) {              //Performing pulling step by the help of second grid (Destination grid)
				if (grid_flag(i,j,0) == 0) {
             		grid_pdf_dest(i,j,0) =grid_pdf(i,j,0);
					grid_pdf_dest(i,j,1) =grid_pdf(i,j-1,1);
					grid_pdf_dest(i,j,2) =grid_pdf(i,j+1,2);
					grid_pdf_dest(i,j,3) =grid_pdf(i+1,j,3);
					grid_pdf_dest(i,j,4) =grid_pdf(i-1,j,4);
					grid_pdf_dest(i,j,5) =grid_pdf(i+1,j-1,5);
					grid_pdf_dest(i,j,6) =grid_pdf(i-1,j-1,6);
					grid_pdf_dest(i,j,7) =grid_pdf(i+1,j+1,7);
					grid_pdf_dest(i,j,8) =grid_pdf(i-1,j+1,8);

				}

			}
		}
	}
//===================================================== VTK Writer ==============================================================
	void write_vtk(const string vtk_file,const Density_Field& grid_density, const Velocity_Field& grid_velocity, const Flag_Field& grid_flag,uint gridsizex, uint gridsizey, uint i){

		string file= vtk_file+to_string(i)+".vtk";
			ofstream ofso(file);
		ofso << "# vtk DataFile Version 3.0" << endl;
		ofso << "SiWiRVisFile" << endl;
		ofso << "ASCII" << endl;
		ofso << "DATASET STRUCTURED_POINTS" << endl;
		ofso << "DIMENSIONS " << gridsizex-2 <<" "<<gridsizey-2<<" "<< 1 <<endl;    
		ofso << "ORIGIN " << 0 <<" "<<0<<" "<<0<<endl;
		ofso << "SPACING "<< 1 <<" "<<1<<" "<<1<<endl;
		ofso << "POINT_DATA "<<(gridsizex-2)*(gridsizey-2)<<endl;                    
		ofso << "SCALARS flags unsigned_int 1"<<endl;
		ofso <<"LOOKUP_TABLE default"<<endl;
		for (uint j =1; j< gridsizey-1; ++j) {                              
			for(uint i = 1; i < gridsizex-1; ++i) {
				ofso << grid_flag(i,j,0) <<endl;
			}
		}	
		ofso << "SCALARS density double 1"<<endl;
		ofso <<"LOOKUP_TABLE default"<<endl;
		for (uint j = 1; j< gridsizey-1; ++j) {                       
			for(uint i = 1; i < gridsizex-1; ++i) {
				ofso << grid_density(i,j,0) <<endl;
		
			}
		}

		ofso << "VECTORS velocity double"<<endl;
		for (uint j = 1; j< gridsizey-1; ++j) {
			for(uint i = 1; i < gridsizex-1; ++i) {
			ofso << grid_velocity(i,j,0) <<" "<<grid_velocity(i,j,1)<<" "<<0<<endl;
			}
		}	
	}

}
class Filereader{
	public:
		void readfile(const std::string filename){
			std::ifstream paramfile(filename);
			if(!paramfile.is_open()){
				std::cout<<"invalid file"<<std::endl;
				exit(1);
			}
			std::string line;
			std::string sparam;
			std::string svalue;
			while(std::getline(paramfile, line)){
				std::istringstream inputline(line);
				inputline>>sparam>>svalue;
				param_map.insert({sparam,svalue});
			}
		}
		template<typename Type>
		Type getparam(const std::string sparam){
			std::istringstream svalue(param_map[sparam]);
			Type value;
			svalue >> value;
			return value;
		}

	private:
		std::map<std::string, std::string> param_map;
};

//===============================================LBM Algorithm============================================================
int main(int argc, char* argv[]) {
	if(argc<2){
		cout<<"No paramter file specified"<<endl;
		return 0;
	}
	string filename = argv[1];
	Filereader reader;
	reader.readfile(filename);
	using namespace lbm;
	uint nx = reader.getparam<uint>("sizex");
	uint ny = reader.getparam<uint>("sizey");
	uint gridsizex = nx + 2;
	uint gridsizey = ny + 2;
	uint timesteps = reader.getparam<uint>("timesteps");
	double uin = reader.getparam<double>("uin");
	uint Re = reader.getparam<uint>("Re");
	double v = uin * ny / Re;
	double t = 3 * v + 0.5;
	uint spherex = reader.getparam<uint>("spherex");
	uint spherey = reader.getparam<uint>("spherey");
	uint diameter = reader.getparam<uint>("diameter");
	string vtk_file = reader.getparam<string>("vtk_file");
	uint vtk_step = reader.getparam<uint>("vtk_step");

	Flag_Field grid_flag(gridsizex, gridsizey);
	PDF_Field grid_pdf(gridsizex, gridsizey);
	PDF_Field grid_pdf_dest(gridsizex, gridsizey);
    PDF_Field grid_pdf_eq(gridsizex,gridsizey);
	Density_Field grid_density(gridsizex,gridsizey);
    Velocity_Field grid_velocity(gridsizex,gridsizey);
	Flag_Initialization(grid_flag, gridsizex, gridsizey, spherex, spherey, diameter);
    PDF_Initialization(grid_pdf, gridsizex, gridsizey, spherex, spherey, diameter);

//================================================== Perfoming Time steps for LBM method ====================================

	for (uint i=0; i<=timesteps; ++i){
		Density(grid_density, grid_pdf,grid_flag, gridsizex, gridsizey);
	    Velocity(grid_velocity, grid_pdf,grid_flag,gridsizex,gridsizey);
		PDF_Equilibrium(grid_pdf_eq,grid_density,grid_velocity,grid_flag,gridsizex,gridsizey);
		Collide (grid_pdf,grid_pdf_eq,grid_flag,t,gridsizex,gridsizey);
		Boundary_handling(grid_pdf,grid_flag,grid_velocity,uin,gridsizex,gridsizey,spherex,spherey,diameter);
		Pulling(grid_pdf_dest,grid_pdf,grid_flag,gridsizex,gridsizey);
		swap(grid_pdf,grid_pdf_dest);
		if (i%vtk_step==0){
			cout<< i <<endl;
			write_vtk(vtk_file,grid_density,grid_velocity,grid_flag,gridsizex,gridsizey,i);
		}
	}

	return 0;
}