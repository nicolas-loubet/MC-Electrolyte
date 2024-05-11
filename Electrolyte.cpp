#include "lib/Configuration.h"
#include <sstream>

const float MOLARITY= 0.001;
const int N_PART= 500;
const int TEMP= 300;
const float SPHERE_R= .5;
Vector* bounds;

Vector* copyVector(Vector* v) {
    return new Vector(v->x, v->y, v->z);
}

void g_r(Configuration** cs, const int N_CONF, const string file_dir) {
    const string CSV_SEP= ",";
    ArrInt bins_pp;
    ArrInt bins_nn;
    ArrInt bins_pn;
    ArrInt bins_np;
    const float DR= .2;
    const float MAX_BINS= int(cs[0]->getBounds()->x/2);
    cout << MAX_BINS/DR << endl;
    bins_pp.size= int(MAX_BINS/DR);
    bins_nn.size= int(MAX_BINS/DR);
    bins_pn.size= int(MAX_BINS/DR);
    bins_np.size= int(MAX_BINS/DR);

    bins_pp.arr= new int[bins_pp.size];
    bins_nn.arr= new int[bins_nn.size];
    bins_pn.arr= new int[bins_pn.size];
    bins_np.arr= new int[bins_np.size];
    for(int i= 0; i < bins_pp.size; i++) {
        bins_pp.arr[i]= 0;
        bins_nn.arr[i]= 0;
        bins_pn.arr[i]= 0;
        bins_np.arr[i]= 0;
    }

    int counter_pp= 0;
    int counter_nn= 0;
    int counter_pn= 0;
    int counter_np= 0;
    for(int i_conf= 0; i_conf < N_CONF; i_conf++) {
        Configuration* c= cs[i_conf];
        for(int i= 1; i <= c->getNParticle(); i++)
            for(int j= 1; j <= c->getNParticle(); j++) {
                if(i == j) continue;
                if(c->getParticle(i)->getCharge() > 0)
                    if(c->getParticle(j)->getCharge() > 0) {
                        bins_pp.arr[getBinPosition(c->getParticle(i)->distanceTo(c->getParticle(j),c->getBounds()),0.,MAX_BINS,bins_pp.size,true)]++;
                        counter_pp++;
                    } else {
                        bins_pn.arr[getBinPosition(c->getParticle(i)->distanceTo(c->getParticle(j),c->getBounds()),0.,MAX_BINS,bins_pn.size,true)]++;
                        counter_pn++;
                    }
                else
                    if(c->getParticle(j)->getCharge() > 0) {
                        bins_np.arr[getBinPosition(c->getParticle(i)->distanceTo(c->getParticle(j),c->getBounds()),0.,MAX_BINS,bins_np.size,true)]++;
                        counter_np++;
                    } else {
                        bins_nn.arr[getBinPosition(c->getParticle(i)->distanceTo(c->getParticle(j),c->getBounds()),0.,MAX_BINS,bins_nn.size,true)]++;
                        counter_nn++;
                    }
            }
    }

    float mean_density= cs[0]->getNParticle()/(2*cs[0]->getBounds()->x*cs[0]->getBounds()->y*cs[0]->getBounds()->z);
    const int N_ION= N_PART/2;
    ofstream f(file_dir);
    f << "r" << CSV_SEP << "gr_++" << CSV_SEP << "gr_--" << CSV_SEP << "gr_+-" << CSV_SEP << "gr_-+" << endl;
    for(int i= 0; i < bins_pp.size; i++) {
        float r= i*DR + DR/2;
        float volume_factor= 4*M_PI*r*r*DR;
        f << r << CSV_SEP << float(bins_pp.arr[i]*N_ION)/(counter_pp*volume_factor*mean_density);
        f << CSV_SEP << float(bins_nn.arr[i]*N_ION)/(counter_nn*volume_factor*mean_density);
        f << CSV_SEP << float(bins_pn.arr[i]*N_ION)/(counter_pn*volume_factor*mean_density);
        f << CSV_SEP << float(bins_np.arr[i]*N_ION)/(counter_np*volume_factor*mean_density) << endl;
    }
    f.flush();
    f.close();

    delete(bins_pp.arr);
    delete(bins_nn.arr);
    delete(bins_pn.arr);
    delete(bins_np.arr);
}

int getNumberConfigurations(const string file_name) {
    ifstream f(file_name);
    string s;
    int counter= 0;
    while(getline(f, s)) {
        int n_atoms= stoi(s);
        for(int i= 0; i <= n_atoms; i++) //comment+atoms
            getline(f,s);
        counter++;
    }
    return counter;
}

Configuration** read_xyz(const string file_name) {
    ifstream f(file_name);
    if(!f) {
        cerr << "File "+file_name+" could not be opened." << endl;
        return nullptr;
    }

    int n_conf= getNumberConfigurations(file_name);
    Configuration** output= new Configuration*[n_conf];
    string type;

    for(int i_conf= 0; i_conf < n_conf; i_conf++) {
        int n_partic;
        string comment;
        string line;
        f >> n_partic;
        getline(f, comment);

        output[i_conf]= new Configuration(n_partic,SPHERE_R,TEMP,copyVector(bounds),78.);

        for(int i_part= 0; i_part < n_partic; i_part++) {
            getline(f, line);
            istringstream iss(line);
            float x,y,z;
            iss >> type >> x >> y >> z;
            Vector* pos= new Vector(x,y,z);
            delete(output[i_conf]->getParticle(i_part+1)->getPosition());
            output[i_conf]->getParticle(i_part+1)->setPosition(pos);
        }

    }

    return output;
}

int main() {
    srand(time(NULL));

    const float L= cbrt(float(N_PART)/(MOLARITY*1.20442816));
    bounds= new Vector(L,L,L);

    const int N_CONF= 20000;
    const float print_every= 100;
    const float escalado= 5;

    remove("conf.xyz");
    remove("eq.xyz");
    remove("prod.xyz");
    remove("gr.csv");

    Configuration* c= new Configuration(N_PART,SPHERE_R,TEMP,copyVector(bounds),78.);
    c->printXYZ("conf.xyz","Random conf");
    float r= c->equilibrate(100000*escalado,1.,5000*escalado);
    c->printXYZ("eq.xyz","Equilibration");
    c->production("prod.xyz",print_every*N_CONF*escalado,print_every*escalado,r);

    delete(c);

    Configuration** cs= read_xyz("prod.xyz");
    g_r(cs,N_CONF,"gr.csv");
    for(int i= 0; i < N_CONF; i++)
        delete(cs[i]);
    delete(cs);

    cout << "Para VMD:" << endl;
    cout << "draw line {0 0 0} {"+to_string(L)+" 0 0}" << endl;
    cout << "draw line {0 0 0} {0 "+to_string(L)+" 0}" << endl;
    cout << "draw line {"+to_string(L)+" 0 0} {"+to_string(L)+" "+to_string(L)+" 0}" << endl;
    cout << "draw line {0 "+to_string(L)+" 0} {"+to_string(L)+" "+to_string(L)+" 0}" << endl;
    cout << "draw line {0 0 "+to_string(L)+"} {"+to_string(L)+" 0 "+to_string(L)+"}" << endl;
    cout << "draw line {0 0 "+to_string(L)+"} {0 "+to_string(L)+" "+to_string(L)+"}" << endl;
    cout << "draw line {"+to_string(L)+" 0 "+to_string(L)+"} {"+to_string(L)+" "+to_string(L)+" "+to_string(L)+"}" << endl;
    cout << "draw line {0 "+to_string(L)+" "+to_string(L)+"} {"+to_string(L)+" "+to_string(L)+" "+to_string(L)+"}" << endl;
    cout << "draw line {0 0 0} {0 0 "+to_string(L)+"}" << endl;
    cout << "draw line {"+to_string(L)+" 0 0} {"+to_string(L)+" 0 "+to_string(L)+"}" << endl;
    cout << "draw line {"+to_string(L)+" "+to_string(L)+" 0} {"+to_string(L)+" "+to_string(L)+" "+to_string(L)+"}" << endl;
    cout << "draw line {0 "+to_string(L)+" 0} {0 "+to_string(L)+" "+to_string(L)+"}" << endl;

    return 0;
}
