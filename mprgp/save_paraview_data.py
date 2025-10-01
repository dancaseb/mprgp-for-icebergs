from paraview.simple import *
import sys

# Get file path from command line argument
if len(sys.argv) > 1:
    vtu_file_path = sys.argv[1]
    output_csv_path = sys.argv[2] if len(sys.argv) > 2 else "output.csv"
else:
    print("No file path provided. Please provide the path to the .vtu file as a command line argument.")
    sys.exit(1)

# Load your data file
data = OpenDataFile(vtu_file_path)

# Apply filters if needed, e.g., Contour, Slice
# contour = Contour(Input=data)
# contour.ContourBy = ['POINTS', 'pressure']

# Update the pipeline to make sure data is processed
Show(data)
RenderAllViews()

# Save data to CSV
SaveData(output_csv_path, proxy=data)
