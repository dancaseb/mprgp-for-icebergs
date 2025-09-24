from paraview.simple import *

# Load your data file
data = OpenDataFile("../square/case_t0001.vtu")

# Apply filters if needed, e.g., Contour, Slice
# contour = Contour(Input=data)
# contour.ContourBy = ['POINTS', 'pressure']

# Update the pipeline to make sure data is processed
Show(data)
RenderAllViews()

# Save data to CSV
SaveData("output.csv", proxy=data)
