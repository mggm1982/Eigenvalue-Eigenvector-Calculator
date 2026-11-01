#include<iostream>
#include<cmath>
using namespace std;
double mat[3][3];

double f(double a,double b,double c,double x){
    return x*x*x+a*x*x+b*x+c;
}

void rowswap(double a[3][3],int row1,int row2){
    for(int i=0;i<3;i++){
        double temp=a[row1][i];
        a[row1][i]=a[row2][i];
        a[row2][i]=temp;
    }
}

void rowreduce(double a[3][3]){
     int row=0;

    for(int col=0; col<3 && row<3; col++)
    {
        int pivot=row;

        for(int i=row+1; i<3; i++)
        {
            if(fabs(a[i][col]) > fabs(a[pivot][col]))
            {
                pivot=i;
            }
        }

        if(fabs(a[pivot][col]) < 0.000001)
        {
            continue;
        }

        rowswap(a,row,pivot);

        double p=a[row][col];

        for(int j=0; j<3; j++)
        {
            a[row][j]=a[row][j]/p;
        }

        for(int i=0; i<3; i++)
        {
            if(i!=row)
            {
                double factor=a[i][col];

                for(int j=0; j<3; j++)
                {
                    a[i][j]=a[i][j]-factor*a[row][j];
                }
            }
        }

        row++;
    }

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            if(fabs(a[i][j])<0.000001)
            {
                a[i][j]=0;
            }
        }
    }
}

void eigenvector(double a[3][3], double k)
{
    double m[3][3];

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            m[i][j]=a[i][j];

            if(i==j)
            {
                m[i][j]=m[i][j]-k;
            }
        }
    }

    rowreduce(m);

    // Find pivot columns
    int pivotcol[3]={-1,-1,-1};
    int rank=0;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            if(fabs(m[i][j])>0.000001)
            {
                pivotcol[rank]=j;
                rank++;
                break;
            }
        }
    }

    // Mark pivot columns
    bool ispivot[3]={false,false,false};

    for(int i=0;i<rank;i++)
    {
        ispivot[pivotcol[i]]=true;
    }

    // Find free columns
    int freecol[3];
    int freecount=0;

    for(int j=0;j<3;j++)
    {
        if(ispivot[j]==false)
        {
            freecol[freecount]=j;
            freecount++;
        }
    }

    cout<<"For eigen value "<<k<<":\n";
    cout<<"Geometric multiplicity = "<<freecount<<'\n';

    // Generate eigenvectors
    for(int f=0;f<freecount;f++)
    {
        double v[3]={0,0,0};

        // Set one free variable to 1
        v[freecol[f]]=1;

        // Calculate pivot variables
        for(int i=0;i<rank;i++)
        {
            int pc=pivotcol[i];

            double sum=0;

            for(int j=0;j<3;j++)
            {
                if(!ispivot[j])
                {
                    sum=sum+m[i][j]*v[j];
                }
            }

            v[pc]=-sum;
        }

        cout<<"Eigenvector "<<f+1<<" = [";

        for(int i=0;i<3;i++)
        {
            cout<<v[i];

            if(i<2)
                cout<<",";
        }

        cout<<"]\n";
    }
}

void cubicroot(double x2,double x1,double x0){
    double lowerbrac=-100;
    for(double i=-99.9;i<=100;i=i+0.1){
        
         if(f(x2,x1,x0,i)*f(x2,x1,x0,lowerbrac) <0){
            double left,right;
            left=lowerbrac;
                right=i;
              double mid;
            for(int j=0;j<100;j++){
                 mid=(left+right)/2;
                if(f(x2,x1,x0,mid)==0){
                    break;
                }
                else if(f(x2,x1,x0,left)*f(x2,x1,x0,mid)<0){
                    right=mid;
                }
                else if(f(x2,x1,x0,right)*f(x2,x1,x0,mid)<0){
                    left=mid;
                }
            }
            cout<<"eigen value is:"<<mid<<'\n';
            eigenvector(mat,mid);
            
        }
        lowerbrac=i;
    }

}
void matprint(double a[3][3],int n){
    cout<<"this is your martix"<<'\n';
    for(int i=0;i<n;i++){
        cout<<"|"<<" ";
        for(int j=0;j<n;j++){
            cout<<a[i][j]<<" ";
        }
        cout<<"|"<<'\n';
    }
}
void eigenvector2x2(double a[3][3],double k){
    double a00;
    double a11;
            a00=a[0][0]-k;
            a11=a[1][1]-k;
            cout<<"for eigen value "<<k<<" the eigen vector is:";
            if (a11==0 && a00==0){
                if(a[0][1]==0){
                    cout<<"[0,1]";
                }
                else if(a[1][0]==0){
                    cout<<"[1,0]";
                }
            }
            else if(a11==0){
                cout<<"["<<-a[0][1]/a00<<",1]\n";
            }
            else if (a00==0){
                cout<<"[1,"<<-a[1][0]/a11<<"]\n";
            }
            else if(a[0][1]==0 && a[1][1]!=0){
                cout<<"[1,"<<-a[1][0]/a11<<"]\n";
            }
            else if(a[1][0]==0){
                cout<<"[1,"<<-a00/a[0][1]<<"]\n";
            }
            
            else{
                cout<<"[1,"<<-a00/a[0][1]<<"]\n";
            }
}
int main(){
    here:
    int order;
    cout<<"type in the order of square matrix(2 or 3):";
    cin >> order;
    if(order==2){
        cout<<"enter the 4 elements of the matrix with space between them: ";
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                cin>>mat[i][j];
            }
        }
        cout<<'\n';
        matprint( mat,2);
        cout<<'\n';
        double y1,y2,trac,descri,deter;
        trac=mat[0][0]+mat[1][1];
        deter=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
        descri=trac*trac-4*deter;
        if(descri<0){
            cout<<"the matrix has imagenary eigen values";
        }
        else{
            y1=(trac+sqrt(descri))/2;
            y2=(trac-sqrt(descri))/2;
            cout<<"eigen values of the matrix are:"<<y1<<" and "<<y2<<'\n';
            eigenvector2x2(mat,y1);
        eigenvector2x2(mat,y2);
        }
        


    }
    else if(order==3){
        cout<<"enter the 9 elements of the matrix with space between them: ";
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                cin>>mat[i][j];
            }
        }
        matprint(mat,3);
        cout<<'\n';
        double trace;
        double adjA;
        double deter;
        deter= mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])- mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])+ mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
        trace=mat[0][0]+mat[1][1]+mat[2][2];
        adjA=(trace*trace-(mat[0][0]*mat[0][0]+ mat[0][1]*mat[1][0]+ mat[0][2]*mat[2][0] + mat[1][0]*mat[0][1] + mat[1][1]*mat[1][1]+ mat[1][2]*mat[2][1] + mat[2][0]*mat[0][2] + mat[2][1]*mat[1][2] + mat[2][2]*mat[2][2]))/2;
        cubicroot(-trace,adjA,-deter);

    }
    else{
        cout<<"error,not a valid input";
        goto here;
    }

}