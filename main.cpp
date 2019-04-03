//Archivo que contiene el procedimiento principal de la implementaci�n,
//describiendo los pasos a seguir para la aplicaci�n del FEM al problema
//de transferencia de calor.

#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"

int main()
{
    //Se preparan dos vectores, uno para contener todas las Ks locales de los elementos de la malla
    //y uno para contener todas las bs locales de dichos elementos
    vector<Matrix> localKs;
    vector<Vector> localbs;

    //Se preparan tambi�n las variables para la K y la b globales, y una para las inc�gnitas de
    //temperatura, que es donde se almacenar� la respuesta.
    Matrix K;
    Vector b;
    Vector T;

    //Se coloca primero una introducci�n con todas las caracter�sticas de la implementaci�n
    cout << "IMPLEMENTACI" << char(224) << "N DEL M" << char(144) << "TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- TRANSFERENCIA DE CALOR\n"
         << "\t- 1 DIMENSI" << char(224) << "N\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n"
         << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    //Se crea un objeto mesh, que contendr� toda la informaci�n de la malla
    mesh m;
    //Se procede a obtener toda la informaci�n de la malla y almacenarla en m
    leerMallayCondiciones(m);

    //Se procede a crear la K local y la b local de cada elemento, almacenando estas
    //estructuras en los vectores localKs y localbs
    crearSistemasLocales(m, localKs, localbs);
    //Descomentar la siguiente l�nea para observar las Ks y bs creadas
    //showKs(localKs); showbs(localbs);

    //Se inicializan con ceros la K global y la b global
    zeroes(K, m.getSize(NODES));
    zeroes(b, m.getSize(NODES));
    //Se procede al proceso de ensamblaje
    ensamblaje(m, localKs, localbs, K, b);
    //Descomentar la siguiente l�nea para observar las estructuras ensambladas
    //showMatrix(K); showVector(b);

    //Se aplican primero las condiciones de contorno de Neumann
    applyNeumann(m, b);
    //Descomentar la siguiente l�nea para observar los cambios en b
    //showVector(b);

    //Luego se aplican las condiciones de contorno de Dirichlet
    applyDirichlet(m, K, b);
    //Descomentar la siguiente l�nea para observar el SEL final luego
    //de los cambios provocados por Dirichlet
    //showMatrix(K); showVector(b);

    //Se prepara con ceros el vector T que contendr� las respuestas
    zeroes(T, b.size());
    //Finalmente, se procede a resolver el SEL
    calculate(K, b, T);

    //Se informa la respuesta:
    cout << "La respuesta es: \n";
    showVector(T);

    return 0;
}
