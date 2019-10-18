/******************************************************************************

SIS dynamics on network
Statistical Mechanics Of Complex Systems
Physics Of Data
Assignment 2
Perin,Polato,Schimmenti

*******************************************************************************/

#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<random>
using namespace std;

//utilities to manipulate C array

float runif(float a= 0, float b = 1)
{
    return a+((b-a)*rand())/RAND_MAX;
}
int runif(int a, int b)
{
    
    return int(runif()*(b-a) + a);
}

float mean(float* x, int length)
{
    float m = 0;
    for(int i = 0; i < length; i++)
    {
        m+=x[i];
    }
    return m/length;
}

float var(float* x, int length)
{
    float m = mean(x, length);
    float v = 0;
    float dx = 0;
    for(int i = 0; i < length; i++)
    {
        dx = x[i]-m;
        v += dx*dx;
    }
    return v/length;
}

struct node
{
    bool s;
    node** sib;
    int deg;
    node()
    {
        s = false;
        sib = NULL;
        deg = 0;
    };
    node(bool state)
    {
        s = state;
        sib = NULL;
        deg = 0;
    };
    node(bool state, int degree)
    {
        s = state;
        sib = new node*[degree];
        deg = degree;
    };
    void init_sibilings(int degree)
    {
        sib = new node*[degree];
        deg = degree;
    }
    void set_sibiling(node* n, int i)
    {
        sib[i] = n;
    };
};


class netlist
{
public:
    node* nodes;
    int N;
    int iCount;
    int ills;
    netlist(int* adjMatrix, int* degrees, int nCount)
    {
        N=nCount;
        nodes = new node[nCount];
        for(int i = 0; i < nCount; i++)
        {
            //cout << "Reading " << i << " node\n";
            int k = degrees[i];
            //cout << "Degree of " << i << " is " << k << "\n";
            nodes[i].init_sibilings(k);
            int h = 0;
            for(int j = 0; j < nCount; j++)
            {
                if(adjMatrix[i*nCount+j] !=0)
                {
                    nodes[i].set_sibiling(nodes+j, h);
                    h++;
                }
            }
            //cout << "Found degree of " << i << " is " << h << "\n---------------------\n";
        }
    }

    void init(int* illInitial, int startIlls)
    {
        for(int i = 0; i < N; i++)
        {
            nodes[i].s = false;
        }
        for(int i = 0; i < startIlls; i++)
        {
            nodes[illInitial[i]].s=true;
        }
        ills = startIlls;
        init_iCount();
    }

    int whereth(bool state, int ith)
    {
        int j = 0;
        int i = 0;
        while(j < ith+1 & i < N)
        {
            if(nodes[i].s==state)
            {
                j++;
            }
            i++;
        }
        if(j < ith)
        {
            return -1;
        }
        else
        {
            return nodes[i-1].s == state ? i-1 : -1;
        }
    }

    int first(bool state)
    {
        for(int i = 0; i < N; i++)
        {
            if(nodes[i].s==state)
            {
                return i;
            }
        }
        cout << "FAIL\n";
        return -1;
    }


    void init_iCount()
    {
        iCount = 0;
        for(int i = 0; i < N; i++)
        {
            bool s_i = nodes[i].s;
            if(s_i)
            {
                int incr = 0;
                for(int j = 0; j < nodes[i].deg; j++)
                {
                    incr += nodes[i].sib[j]->s ? 0 : 1;
                }
                iCount+=incr;
            }
        }
    }
    void update(int i)
    {
        bool s_i = nodes[i].s;
        int incr = s_i ? -1 : 1;
        ills += -incr;
        for(int j = 0; j < nodes[i].deg; j++)
        {
            iCount += nodes[i].sib[j]->s ? incr : -incr;
        }
    }

    float simulate_sis(float lmbd, float mu, int iters, int maxInfected)
    {
        int* indices = new int[maxInfected];
        for(int i = 0; i < maxInfected; i++)
        {
            indices[i] = runif(0, N);
        }
        init(indices, maxInfected);
        //cout << "Simulation starting with " << ills << " infected and susceptible edges " << iCount << "\n";
        for(int t = 0; t < iters; t++)
        {
            if(ills == 0 && iCount>0)
            {
                cout << "IMPOSSIBLE\n";
                return 0;   
            }
            if(ills < 0 || iCount < 0)
            {
                cout << "FATAL\n";
                return 0;
            }
            
            if(ills == 0 || iCount==0)
            {
                //cout << "Finished due to ill death...\n";
                break;
            }
            
            float p_rec = ills*mu/(ills*mu+iCount*lmbd);

            bool recover = runif()<p_rec;  
            int rndPos = runif(0, recover ? ills : N-ills);
            //int i = whereth(recover, rndPos);
            int i = first(recover);
            if(i >= 0)
            {
                //cout << "Ills: " << ills << "\n";
                nodes[i].s = !recover;
                update(i);
            }
            
            
        }
        float inf = 0;
        for(int i = 0; i < N; i++)
        {
            inf += nodes[i].s ? 1 : 0;
        }
        return inf/N;
    }

};



int main(int argc, char *argv[]) //fname, nCount, lmbdFrom, lmbdTo, lmbdCount, mu, iters, nSamples, maxInfected
{
    try
    {
        
    
    char* fName = argv[1];
    cout << "Filename: " << fName << "\n";
    int nCount = stoi(argv[2]);
    cout << "Node Count: " << nCount << "\n";
    float lmbdFrom = stof(argv[3]);
    cout << "Lambda from: " << lmbdFrom << "\n";
    float lmbdTo = stof(argv[4]);
    cout << "Lambda to: " << lmbdTo << "\n";
    int lmbdCount = stoi(argv[5]);
    cout << "Lambda count: " << lmbdCount << "\n";
    float mu = stof(argv[6]);
    cout << "Mu: " << mu << "\n";
    int iters = stoi(argv[7]);
    cout << "Iters: " << iters << "\n";
    int nSamples = stoi(argv[8]);
    cout << "nSamples: " << nSamples << "\n";
    int maxInfected = stoi(argv[9]);
    cout << "maxInfected: " << maxInfected << "\n";
    char* outName = argv[10];
    cout << "Output Filename: " << outName << "\n";
    ifstream inFile(fName);
    int i,j;
    int nVerb = 0;
    int* adjMatrix = new int[nCount*nCount];
    int* degrees = new int[nCount];
    for(int i = 0; i < nCount; i++)
    {
        degrees[i]=0;
    }
    while(inFile >> i >> j)
    {
        adjMatrix[i*nCount+j] = 1;
        adjMatrix[j*nCount+i] = 1;
        //cout << "(" << i << "," << j << ")\n";
    }
    for(int i = 0; i < nCount; i++)
    {
        for(int j= 0; j < nCount; j++)
        {
            degrees[i]+=adjMatrix[i*nCount+j];
        }
    }
    inFile.close();
    cout << "Creating the network...\n";
    netlist net(adjMatrix, degrees, nCount);
    cout << "Network created...\n";
    float l = lmbdFrom;
    float step = (lmbdTo-lmbdFrom)/lmbdCount;
    float* infs = new float[lmbdCount];
    ofstream outFile(outName);
    for(int k = 0; k < lmbdCount; k++)
    {
        cout << "Starting...\n";
        
        try
        {
            outFile << l << "\t";
            float mInf = 0;
            for(int n = 0; n < nSamples; n++)
            {
                float inf =net.simulate_sis(l, mu, iters, maxInfected);
                if(n==nSamples-1)
                {
                    outFile << inf << "\n";
                }
                else
                {
                    outFile << inf << "\t";
                }
                mInf +=inf;
                if(nVerb != 0 && n%nVerb==0)
                {
                    cout << "Lambda: " << l << " infected: " << inf << " and sample " << n << "\n";
                }
            }
            mInf/=nSamples;
            infs[k] = mInf;
            cout << "Lambda: " << l << " Infected: " << mInf << "\n";
            l+=step;
        }
        catch (...)
        {
            outFile.close();
            cout << "Error\n";
        }
        cout << "Finishing...\n--------------------------\n";
    }
    outFile.close();

    
    }
    catch(...)
    {
        cout << "Call with the arguments: fname, nCount, lmbdFrom, lmbdTo, lmbdCount, mu, iters, nSamples, maxInfected, outFile\n";
    }
}
