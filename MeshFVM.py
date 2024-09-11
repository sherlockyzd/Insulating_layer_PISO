import numpy as np

class MeshStructure:
    def __init__(self, dimension, dims, cellsize, cellcenters, facecenters, corners, edges):
        self.dimension = dimension
        self.dims = dims
        self.cellsize = cellsize
        self.cellcenters = cellcenters
        self.facecenters = facecenters
        self.corners = corners
        self.edges = edges

class MeshGenerator:
    def create_mesh_2d(self, Nx, Ny, Width, Height):
        """
        Create a uniform 2D mesh.
        
        Nx: Number of cells in the x direction
        Ny: Number of cells in the y direction
        Width: Domain length in the x direction
        Height: Domain length in the y direction
        """
        dx = Width / Nx
        dy = Height / Ny

        # Cell centers
        CellLocation_x = np.linspace(dx / 2, Width - dx / 2, Nx)+dx
        CellLocation_y = np.linspace(dy / 2, Height - dy / 2, Ny)

        # Face locations
        FaceLocation_x = np.linspace(0, Width, Nx + 1)+dx
        FaceLocation_y = np.linspace(0, Height, Ny + 1)

        # Cell sizes
        CellSize_x = np.ones(Nx ) * dx
        CellSize_y = np.ones(Ny ) * dy

        # Generate grid numbering
        G = np.reshape(np.arange(1, (Nx ) * (Ny ) + 1), (Nx , Ny ))

        corners = G[[0, -1], [0, -1]].flatten()

        return MeshStructure(
            dimension=2,
            dims=[Nx, Ny],
            cellsize={'x': CellSize_x, 'y': CellSize_y, 'z': [0.0]},
            cellcenters={'x': CellLocation_x, 'y': CellLocation_y, 'z': [0.0]},
            facecenters={'x': FaceLocation_x, 'y': FaceLocation_y, 'z': [0.0]},
            corners=corners,
            edges=[1]
        )

    def create_radial_mesh_2d(self, Nr, Ntetta, R, Tetta):
        """
        Create a radial 2D mesh (for circular or cylindrical domains).
        
        Nr: Number of cells in the radial direction
        Ntetta: Number of cells in the angular direction
        R: Radius of the domain
        Tetta: Angular span of the domain (up to 2*pi)
        """
        if Tetta > 2 * np.pi:
            Tetta = 2 * np.pi
            print("Warning: Tetta was greater than 2*pi, and it has been scaled down to 2*pi.")

        return self.create_mesh_2d(Nr, Ntetta, R, Tetta)

    def create_nonuniform_mesh_2d(self, facelocationX, facelocationY):
        """
        Create a non-uniform 2D mesh based on specified face locations.
        
        facelocationX: Face locations in the x direction
        facelocationY: Face locations in the y direction
        """
        facelocationX = np.array(facelocationX).flatten()
        facelocationY = np.array(facelocationY).flatten()

        Nx = len(facelocationX) - 1
        Ny = len(facelocationY) - 1

        # Cell centers and sizes
        CellLocation_x = 0.5 * (facelocationX[1:] + facelocationX[:-1])
        CellLocation_y = 0.5 * (facelocationY[1:] + facelocationY[:-1])

        CellSize_x = np.diff(facelocationX)
        CellSize_y = np.diff(facelocationY)

        # Generate grid numbering
        G = np.reshape(np.arange(1, (Nx + 2) * (Ny + 2) + 1), (Nx + 2, Ny + 2))

        corners = G[[0, -1], [0, -1]].flatten()

        return MeshStructure(
            dimension=2,
            dims=[Nx, Ny],
            cellsize={'x': CellSize_x, 'y': CellSize_y, 'z': [0.0]},
            cellcenters={'x': CellLocation_x, 'y': CellLocation_y, 'z': [0.0]},
            facecenters={'x': facelocationX, 'y': facelocationY, 'z': [0.0]},
            corners=corners,
            edges=[1]
        )
