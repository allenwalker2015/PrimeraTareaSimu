//La funci�n recibe:
//- Una matriz
//La funci�n muestra en pantalla el contenida de la matrix, fila por fila
void showMatrix(Matrix K)
{
    for (int i = 0; i < K.at(0).size(); i++)
    {
        cout << "[\t";
        for (int j = 0; j < K.size(); j++)
        {
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}

//La funci�n recibe un arreglo de matrices, y las muestra
//en pantalla una por una
void showKs(vector<Matrix> Ks)
{
    for (int i = 0; i < Ks.size(); i++)
    {
        cout << "K del elemento " << i + 1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}

//La funci�n recibe:
//- Un vector
//La funci�n asume que recibe un vector columna y muestra su contenido
//en pantalla en una sola fila
void showVector(Vector b)
{
    cout << "[\t";
    for (int i = 0; i < b.size(); i++)
    {
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

//La funci�n recibe un arreglo de vectores, y los muestra
//en pantalla uno por uno
void showbs(vector<Vector> bs)
{
    for (int i = 0; i < bs.size(); i++)
    {
        cout << "b del elemento " << i + 1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

//La funci�n recibe:
//- Un elemento
//- El objeto mesh
//La funci�n construye la matrix local K para el elemento
//especificado de acuerdo a la formulaci�n del problema
Matrix createLocalK(int element, mesh &m)
{
    //Se prepara la matriz y sus dos filas (se sabe que es una matriz 2x2)
    Matrix K;
    Vector row1, row2;

    //De acuerdo a la formulaci�n, la matriz local K tiene la forma:
    //          (k/l)*[ 1 -1 ; -1 1 ]

    float l = m.getNodes()[element + 1].getX() -  m.getNodes()[element].getX();
    cout << "El elemento Xn es:" <<m.getNodes()[element].getX();
    cout << "El elemento Xn+1 es:" <<m.getNodes()[element+1].getX();
    cout << "Local K => La L es:" << l << endl;
    //Se extraen del objeto mesh los valores de k y l
    float k = m.getParameter(THERMAL_CONDUCTIVITY);
    //Se crean las filas
    row1.push_back(k / l);
    row1.push_back(-k / l);
    row2.push_back(-k / l);
    row2.push_back(k / l);
    //Se insertan las filas en la matriz
    K.push_back(row1);
    K.push_back(row2);

    return K;
}

//La funci�n recibe:
//- Un elemento
//- El objeto mesh
//La funci�n construye el vector local b para el elemento
//especificado de acuerdo a la formulaci�n del problema
Vector createLocalb(int element, mesh &m)
{
    //Se prepara el vector b (se sabe que ser� un vector 2x1)
    Vector b;

    //Se sabe que el vector local b tendr� la forma:
    //          (Q*l/2)*[ 1 ; 1 ]

    //Se extraen del objeto mesh los valores de Q y l
    float l = m.getNodes()[element + 1].getX() -  m.getNodes()[element].getX();
    //cout<<"Local B => La L es:"<<l;

    float Q = m.getParameter(HEAT_SOURCE);
    //Se insertan los datos en el vector
    b.push_back(Q * l / 2);
    b.push_back(Q * l / 2);

    return b;
}

//La funci�n recibe:
//- El objeto mesh
//- Un arreglo de matrices
//- Un arreglo de vectores
//La funci�n construye una K y una b locales para cada elemento de la malla,
//y los almacena en su arreglo correspondiente
void crearSistemasLocales(mesh &m, vector<Matrix> &localKs, vector<Vector> &localbs)
{
    //Se recorren los elementos
    for (int i = 0; i < m.getSize(ELEMENTS); i++)
    {
        //Por cada elemento, se crea su K y su b
        localKs.push_back(createLocalK(i, m));
        localbs.push_back(createLocalb(i, m));
    }
}

//La funci�n recibe:
//- El elemento actual
//- La matriz local K
//- La matriz global K
//La funci�n inserta la K local en la K global de acuerdo a los nodos
//del elemento
void assemblyK(element e, Matrix localK, Matrix &K)
{
    //Se determinan los nodos del elemento actual como los �ndices de la K global
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    //Se utilizan los �ndices para definir las celdas de la submatriz
    //a la que se agrega la matriz local del elemento actual
    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);
}

//La funci�n recibe:
//- El elemento actual
//- El vector local b
//- El vector global b
//La funci�n inserta la b local en la b global de acuerdo a los nodos
//del elemento
void assemblyb(element e, Vector localb, Vector &b)
{
    //Se determinan los nodos del elemento actual como los �ndices de la b global
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    //Se utilizan los �ndices para definir las celdas del subvector
    //al que se agrega el vector local del elemento actual
    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
}

//La funci�n recibe:
//- El objeto mesh
//- El arreglo de Ks locales
//- El arreglo de bs locales
//- La matriz K global
//- El vector b global
//La funci�n se encarga de ensamblar adecuadamente todos los sistemas locales en
//la K y la b globales
void ensamblaje(mesh &m, vector<Matrix> &localKs, vector<Vector> &localbs, Matrix &K, Vector &b)
{
    //Se recorren todos los elementos de la malla, uno por uno
    for (int i = 0; i < m.getSize(ELEMENTS); i++)
    {
        //Se extrae del objeto mesh el elemento actual
        element e = m.getElement(i);
        //Se ensamblan la K y la b del elemento actual en las variables globales
        assemblyK(e, localKs.at(i), K);
        assemblyb(e, localbs.at(i), b);
    }
}

//La funci�n recibe:
//- El objeto mesh
//- El vector b global
//La funci�n aplica en la b global las condiciones de Neumann en las
//posiciones que correspondan
void applyNeumann(mesh &m, Vector &b)
{
    //Se recorren las condiciones de Neumann, una por una
    for (int i = 0; i < m.getSize(NEUMANN); i++)
    {
        //Se extrae la condici�n de Neumann actual
        condition c = m.getCondition(i, NEUMANN);
        //En la posici�n de b indicada por el nodo de la condici�n,
        //se agrega el valor indicado por la condici�n
        b.at(c.getNode1() - 1) += c.getValue();
    }
}

//La funci�n recibe:
//- El objeto mesh
//- La matriz K global
//- El vector b global
//La funci�n aplica en la K y b globales las condiciones de Dirichlet, eliminando
//las filas correspondientes, y enviando desde el lado izquierdo del SEL al lado
//derecho los valores de las columnas correspondientes
void applyDirichlet(mesh &m, Matrix &K, Vector &b)
{
    //Se recorren las condiciones de Dirichlet, una por una
    for (int i = 0; i < m.getSize(DIRICHLET); i++)
    {
        //Se extrae la condici�n de Dirichlet actual
        condition c = m.getCondition(i, DIRICHLET);
        //Se establece el nodo de la condici�n como el �ndice
        //para K y b globales donde habr� modificaciones
        int index = c.getNode1() - 1;

        //Se elimina la fila correspondiente al nodo de la condici�n
        K.erase(K.begin() + index); //Se usa un iterator a la posici�n inicial, y se
        b.erase(b.begin() + index); //le agrega la posici�n de inter�s

        //Se recorren las filas restantes, una por una, de modo que
        //el dato correspondiente en cada fila a la columna del nodo de la
        //condici�n, se multiplique por el valor de Dirichlet, y se env�e al
        //lado derecho del SEL con su signo cambiado
        for (int row = 0; row < K.size(); row++)
        {
            //Se extrae el valor ubicado en la columna correspondiente
            //al nodo de la condici�n
            float cell = K.at(row).at(index);
            //Se elimina la columna correspondiente
            //al nodo de la condici�n
            K.at(row).erase(K.at(row).begin() + index);
            //El valor extra�do se multiplica por el valor de la condici�n,
            //se le cambia el signo, y se agrega al lado derecho del SEL
            b.at(row) += -1 * c.getValue() * cell;
        }
    }
}

//La funci�n recibe:
//- La matriz K global
//- El vector b global
//- El vector T que contendr� los valores de las inc�gnitas
//La funci�n se encarga de resolver el SEL del problema
void calculate(Matrix &K, Vector &b, Vector &T)
{
    //Se utiliza lo siguiente:
    //      K*T = b
    // (K^-1)*K*T = (K^-1)*b
    //     I*T = (K^-1)*b
    //      T = (K^-1)*b

    //Se prepara la inversa de K
    Matrix Kinv;
    //Se calcula la inversa de K
    inverseMatrix(K, Kinv);
    //Se multiplica la inversa de K por b, y el resultado se almacena en T
    productMatrixVector(Kinv, b, T);
}
