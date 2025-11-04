#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

ILOSTLBEGIN

struct Edge 
{
    int i;
    int j;
    double g;
    double p;
};

static void usage(const char *progname);

int main(int argc, char **argv) 
{
    IloEnv env;
    try 
    {
        if(argc != 2) 
        {
            usage(argv[0]);
            throw(-1);
        }

        std::ifstream file(argv[1]);
        if(!file)
        {
            std::cerr << "Erro: Não foi possível abrir o arquivo " << argv[1] << std::endl;
            throw(-1);
        }

        int ni, nj, NL;
        double c, Q;
        std::vector<Edge> edges;

        file >> ni >> nj >> c >> Q >> NL;

        for(int k = 0; k < NL; k++) 
        {
            int i, j;
            double g, p;
            file >> i >> j >> g >> p;
            edges.push_back({i - 1, j - 1, g, p});
        }

        file.close();

        IloModel model(env);
        
        IloNumVarArray y(env, ni, 0, 1, ILOBOOL);
        for(int i = 0; i < ni; i++) 
        {
            std::stringstream name;
            name << "y(" << i + 1 << ")";
            y[i].setName(name.str().c_str());
        }

        IloNumVarArray x(env, NL, 0, 1, ILOBOOL);
        for(int k = 0; k < NL; k++) 
        {
            std::stringstream name;
            name << "x(" << edges[k].i+1 << "," << edges[k].j+1 << ")";
            x[k].setName(name.str().c_str());
        }

        IloExpr objExpr(env);
        for(int i = 0; i < ni; i++) 
        {
            objExpr += c * y[i];
        }
        for(int k = 0; k < NL; k++)
        {
            objExpr += edges[k].g * x[k];
        }
        model.add(IloMinimize(env, objExpr));
        objExpr.end();

        for(int j = 0; j < nj; ++j) 
        {
            IloExpr r1Expr(env);
            for(int k = 0; k < NL; k++)
            {
                if(edges[k].j == j) r1Expr += x[k];
            }
            model.add(r1Expr == 1);
            r1Expr.end();
        }

        for(int k = 0; k < NL; k++)
        {
            int i = edges[k].i;
            model.add(x[k] <= y[i]);
        }

        for(int i = 0; i < ni; i++) 
        {
            IloExpr r3Expr(env);
            for(int k = 0; k < NL; k++) 
            {
                if(edges[k].i == i) r3Expr += edges[k].p * x[k];
            }
            model.add(r3Expr <= Q * y[i]);
            r3Expr.end();
        }

        IloCplex cplex(model);

        cplex.solve();

        env.out() << "Status da Solucao = " << cplex.getStatus() << std::endl;
        env.out() << "Valor da Solucao (Custo Minimo) = " << cplex.getObjValue() << std::endl; //

        IloNumArray y_vals(env);
        cplex.getValues(y_vals, y);
        env.out() << "Facilidades Construidas:" << std::endl;
        for(int i = 0; i < ni; i++) 
        {
            if(y_vals[i] > 0.9)
            { 
                env.out() << "  Facilidade " << i + 1 << std::endl;
            }
        }

        IloNumArray x_vals(env);
        cplex.getValues(x_vals, x);
        env.out() << "Atendimentos:" << std::endl;
        for(int k = 0; k < NL; k++) 
        {
            if(x_vals[k] > 0.9) 
            {
                env.out() << "  Facilidade " << edges[k].i + 1 << " -> Cliente " << edges[k].j + 1 << std::endl;
            }
        }
        
        y_vals.end();
        x_vals.end();

    } 
    catch(IloException& e) 
    {
        std::cerr << "Erro no CPLEX: " << e << std::endl;
    } 
    catch(...) 
    {
        std::cerr << "Erro desconhecido" << std::endl;
    }

    env.end();
    return 0;
}

static void usage(const char *progname) 
{
    std::cerr << "Uso: " << progname << " <arquivo_da_instancia.txt>" << std::endl;
    std::cerr << " Ex: " << progname << " TPI_F_0.txt" << std::endl;
}