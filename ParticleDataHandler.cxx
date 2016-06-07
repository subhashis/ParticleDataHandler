#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtkDelaunay3D.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vector>
#include <string>

#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>

using namespace std;
/*Density field creator for particle data*/
void createDensityField(int dataID){
    string number;
    if(dataID < 10)
        number = "00" + std::to_string(dataID);
    else if(dataID < 100)
        number = "0" + std::to_string(dataID);
    else
        number = std::to_string(dataID);

    std::string filename = "/media/WD3/SciVisContest16/smoothinglength_0.20/run30/" + number +".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkDataSet *readData = reader->GetOutputAsDataSet();
    /*std:cout << "subhashis\n";
    std::cout << "no. of cells=" << readData->GetNumberOfCells() << std::endl;
    std::cout << "no. of points=" << readData->GetNumberOfPoints() << std::endl;
    std::cout << "no. of arrays=" << readData->GetPointData()->GetNumberOfArrays() << std::endl;
    std::cout << "name=" << readData->GetPointData()->GetArray(0)->GetName() << std::endl;
    std::cout << "maxId=" << readData->GetPointData()->GetArray(0)->GetMaxId() << std::endl;
    std::cout << "name=" << readData->GetPointData()->GetArray(1)->GetName() << std::endl;
    std::cout << "maxId=" << readData->GetPointData()->GetArray(1)->GetMaxId() << std::endl;*/

    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(readData->GetPointData()->GetArray("concentration"));

    float maxC = -10000.0;
    float minC = 10000.0;

    int nop = readData->GetNumberOfPoints();

    for(int i = 0; i < nop; i++)
    {
          float value;
          value = array->GetValue(i);
          if(value > maxC) maxC=value;
          if(value < minC) minC=value;
          //std::cout << i << ": " << value << std::endl;
    }
    std::cout << "max = " << maxC << " min = " << minC << std::endl;

    //double a[3];
    //readData->GetPoint(1661832,a);
    //std::cout<< a[0] << " " << a[1] << " " << a[2] << std::endl;

    //find the max-min of x,y,z.
    double xmax = -1000.0, xmin = 1000.0;
    double ymax = -1000.0, ymin = 1000.0;
    double zmax = -1000.0, zmin = 1000.0;
    for(int i=0; i<nop; i++)
    {
        double a[3];
        readData->GetPoint(i,a);

        if(xmax < a[0]) xmax = a[0];
        if(xmin > a[0]) xmin = a[0];

        if(ymax < a[1]) ymax = a[1];
        if(ymin > a[1]) ymin = a[1];

        if(zmax < a[2]) zmax = a[2];
        if(zmin > a[2]) zmin = a[2];
    }

    /*cout << "bounding box limit\n";
    cout << "x:[" << xmin << " - " << xmax << "]" << endl;
    cout << "y:[" << ymin << " - " << ymax << "]" << endl;
    cout << "z:[" << zmin << " - " << zmax << "]" << endl;*/

    double xLen = xmax - xmin;
    double yLen = ymax - ymin;
    double zLen = zmax - zmin;

    int gridResX = 256;
    int gridResY = 256;
    int gridResZ = 256;

    FILE *fOut;
    double *outData;
    string f1 = "/media/WD3/SciVisContest16/smoothinglength_0.20/run30/density" + number + "_x256.raw";
    //fOut = fopen(f1.c_str(),"wb");
    outData = new double[gridResX * gridResY * gridResZ];

    int perGridLenX = xLen / gridResX;
    int perGridLenY = yLen / gridResY;
    int perGridLenZ = zLen / gridResZ;


    std::vector< std::vector< std::vector< std::vector< int > > > > pointsInGrid;
    pointsInGrid.resize(gridResX);
    for (int j = 0; j < gridResX; j++){
        pointsInGrid[j].resize(gridResY);
        for (int k = 0; k < gridResY; k++){
            pointsInGrid[j][k].resize(gridResZ);
        }
    }

    for (int i = 0; i < nop; i++){
        double a[3];
        readData->GetPoint(i,a);

        int xLoc = ((a[0] - xmin) * (gridResX-1)) / xLen;
        int yLoc = ((a[1] - ymin) * (gridResY-1)) / yLen;
        int zLoc = ((a[2] - zmin) * (gridResZ-1)) / zLen;

        pointsInGrid[xLoc][yLoc][zLoc].push_back(i);

    }

    //Now we will do the reverse wighed interpolation
    //std::cout << "[INFO] Started creating the actual grid.\n";

    for (int i = 0; i < gridResX; i++){
        for (int j = 0; j < gridResY; j++){
            for (int k = 0; k < gridResZ; k++){

                //There can be max of 8 cells and min of 6 cells
                std::vector<int> allPoints;

                std::vector<int> grid1 = pointsInGrid[i][j][k];
                allPoints.resize(grid1.size());
                std::copy(grid1.begin(), grid1.end(), allPoints.begin());

                if (j - 1 >= 0) {
                    std::vector<int> grid2 = pointsInGrid[i][j - 1][k];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid2.size();
                    allPoints.resize(newSize);
                    std::copy(grid2.begin(), grid2.end(), allPoints.begin() + curSize);
                }

                if (k - 1 >= 0){
                    std::vector<int> grid3 = pointsInGrid[i][j][k - 1];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid3.size();
                    allPoints.resize(newSize);
                    std::copy(grid3.begin(), grid3.end(), allPoints.begin() + curSize);
                }

                if ((j - 1 >= 0) && (k - 1 >= 0)) {
                    std::vector<int> grid4 = pointsInGrid[i][j - 1][k - 1];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid4.size();
                    allPoints.resize(newSize);
                    std::copy(grid4.begin(), grid4.end(), allPoints.begin() + curSize);
                }

                if (i - 1 >= 0) {
                    std::vector<int> grid5 = pointsInGrid[i - 1][j][k];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid5.size();
                    allPoints.resize(newSize);
                    std::copy(grid5.begin(), grid5.end(), allPoints.begin() + curSize);
                }

                if ((i - 1 >= 0) && (j - 1 >= 0)) {
                    std::vector<int> grid6 = pointsInGrid[i - 1][j - 1][k];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid6.size();
                    allPoints.resize(newSize);
                    std::copy(grid6.begin(), grid6.end(), allPoints.begin() + curSize);
                }

                if ((i - 1 >= 0) && (k - 1 >= 0)){
                    std::vector<int> grid7 = pointsInGrid[i - 1][j][k - 1];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid7.size();
                    allPoints.resize(newSize);
                    std::copy(grid7.begin(), grid7.end(), allPoints.begin() + curSize);
                }

                if ((i - 1 >= 0) && (j - 1 >= 0) && (k - 1 >= 0)) {
                    std::vector<int> grid8 = pointsInGrid[i - 1][j - 1][k - 1];
                    int curSize = allPoints.size();
                    int newSize = curSize + grid8.size();
                    allPoints.resize(newSize);
                    std::copy(grid8.begin(), grid8.end(), allPoints.begin() + curSize);
                }

                double gridX = xmin + (i*perGridLenX);
                double gridY = ymin + (i*perGridLenY);
                double gridZ = zmin + (i*perGridLenZ);

                int totalSize = allPoints.size();
                if (totalSize > nop) {
                    std::cout << "That happened.\n";
                    outData[k * gridResX * gridResY + j * gridResX + i] = 0.0;
                }
                else {
                    double val = 0.0;
                    double sumOfWeights = 0.0;
                    for (std::vector<int>::iterator it = allPoints.begin(); it != allPoints.end(); it++){
                        int idx = *it;

                        double a[3];
                        readData->GetPoint(idx,a);
                        double alpha = array->GetValue(idx);
                        //double alpha = inData[4 * idx + 3];
                        double distSq = ((a[0] - gridX)*(a[0] - gridX)) + ((a[1] - gridY)*(a[1] - gridY)) + ((a[2] - gridZ)*(a[2] - gridZ));
                        double dist = sqrt(distSq);
                        double weight = 1.0 / dist;
                        val += (alpha * weight);
                        sumOfWeights += weight;
                        //std::cout << idx << "  " << weight << "\n";
                    }
                    //getchar();
                    if (totalSize != 0)
                        outData[k * gridResX * gridResY + j * gridResX + i] = val/sumOfWeights;
                    else
                        outData[k * gridResX * gridResY + j * gridResX + i] = 0.0;
                }


            }
        }

    }
    std::cout << number <<  "\nDONE\n";
    //fwrite (outData , sizeof(double), gridResX*gridResY*gridResZ, fOut);
    //fclose(fOut);



    //write to vti file
    vtkSmartPointer<vtkImageData> imageData =  vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(gridResX,gridResY,gridResZ);
    #if VTK_MAJOR_VERSION <= 5
      imageData->SetNumberOfScalarComponents(1);
      imageData->SetScalarTypeToDouble();
    #else
      imageData->AllocateScalars(VTK_DOUBLE, 1);
    #endif
      int* dims = imageData->GetDimensions();

      // Fill every entry of the image data with "2.0"
      for (int z = 0; z < dims[2]; z++)
        {
        for (int y = 0; y < dims[1]; y++)
          {
          for (int x = 0; x < dims[0]; x++)
            {
            double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
            pixel[0] = outData[z * gridResX * gridResY + y * gridResX + x];
            }
          }
        }

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
      string f2 = "/media/WD3/SciVisContest16/smoothinglength_0.20/run30/density" + number + "_x256.vti";
      writer->SetFileName(f2.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInputConnection(imageData->GetProducerPort());
    #else
      writer->SetInputData(imageData);
    #endif
      writer->Write();

}

int main ( int argc, char *argv[] )
{


    //read all the data from the file
    /*for(int i=10; i<60;i++)
        createDensityField(i);*/



    return EXIT_SUCCESS;
}
