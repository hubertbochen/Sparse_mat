class Mesh2D
{
    public:
    
    int nx;
    int ny;
    int xlim;
    int ylim;
    struct Dense pt_type;// non_exist: 0 , boundary:1 , inner :2
    int* xlist;// xlist(lin) to find x index
    int* ylist;// ylist(lin) to find y index
    int* linlist;// linlist(i,j) to find the lin index
    
    
    
    struct SpMat stiff;    
    
    void build_list(){// input: pt_type
    
    };
    
    void build_stiff()
    {
    
    };
    
}