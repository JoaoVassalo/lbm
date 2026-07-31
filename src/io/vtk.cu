#include "vtk.cuh"

__host__ void writeVTI(size_t t, const std::string &outDir, varType *mom_host)
{
    namespace fs = std::filesystem;

    fs::path vtkDir = fs::path(outDir) / "vtk";
    fs::create_directories(vtkDir);

    std::ostringstream filename;
    filename << "output" << std::setw(6) << std::setfill('0') << t << ".vti";

    fs::path filepath = vtkDir / filename.str();

    std::ofstream file(filepath.string());
    if (!file.is_open())
    {
        std::cerr << "Could not open VTI file for writing: " << filepath.string() << "\n";
        return;
    }

    const int nx = Geometry::NX;
    const int ny = Geometry::NY;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <ImageData WholeExtent=\"0 " << (nx - 1)
         << " 0 " << (ny - 1)
         << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    file << "    <Piece Extent=\"0 " << (nx - 1)
         << " 0 " << (ny - 1)
         << " 0 0\">\n";

    file << "      <PointData>\n";

    // rho
    file << "        <DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const varType rho = mom_host[momIdx<momId::rho>(x, y)]; // rho
            file << "          " << static_cast<varType>(rho) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // ux
    file << "        <DataArray type=\"Float32\" Name=\"ux\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const varType ux = mom_host[momIdx<momId::ux>(x, y)]; // ux
            file << "          " << static_cast<varType>(ux) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // uy
    file << "        <DataArray type=\"Float32\" Name=\"uy\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const varType uy = mom_host[momIdx<momId::uy>(x, y)]; // uy
            file << "          " << static_cast<varType>(uy) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // velocity
    file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const varType ux = mom_host[momIdx<momId::ux>(x, y)]; // ux
            const varType uy = mom_host[momIdx<momId::uy>(x, y)]; // uy

            file << "          "
                 << static_cast<varType>(ux) << " "
                 << static_cast<varType>(uy) << " "
                 << 0.0f << "\n";
        }
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
    file << "      <CellData>\n";
    file << "      </CellData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";

    file.close();
}