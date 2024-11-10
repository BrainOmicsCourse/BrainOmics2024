import plotly.express as px
from plotly.colors import sample_colorscale


# Create SankeyDF dataframe (same as your example)



def plotSankey(adata, covs=[]):

	SankeyDF = adata.obs[covs]

	colorMap = {i:adata.uns[covs[0]+"_colors"][adata.obs[covs[0]].cat.categories.tolist().index(i)] for i in adata.obs[covs[0]].unique().tolist()}

	SankeyDF["color"] = SankeyDF[covs[0]].replace(colorMap)


	# Create the parallel categories plot without the colorbar
	fig = px.parallel_categories(SankeyDF,
								dimensions=covs,
								color="color",
								color_continuous_scale=px.colors.sequential.RdBu,
								width=1200, height=700)  # Increased width

	# Remove the colorbar by setting showscale=False
	fig.update_layout(coloraxis_showscale=False)

	# Adjust margins to increase plot area and make room for the legend
	fig.update_layout(
		margin=dict(l=100, r=300, t=50, b=50)  # Increased right margin to 200
	)

	# Get exact color values from the continuous scale
	color_scale = px.colors.sequential.RdBu



	# Add mock legend (colored rectangles with annotations)
	x_start = 1.06  # Move the legend closer to the plot (from 1.02)
	y_start = 0.9   # Start position for the legend
	rect_height = 0.05  # Height of each rectangle in the mock legend

	for i, (label, color) in enumerate(colorMap.items()):
		y_pos = y_start - i * rect_height  # Adjust the y position for each item in the legend

		# Add colored rectangle (the mock color legend)
		fig.add_shape(type="rect",
					xref="paper", yref="paper",
					x0=x_start, y0=y_pos, x1=x_start + 0.03, y1=y_pos + rect_height,
					line=dict(color='black'), fillcolor=color)  # Added black outline for visibility

		# Add corresponding label right next to the rectangle
		fig.add_annotation(x=x_start + 0.035, y=y_pos + rect_height / 2,
						xref="paper", yref="paper",
						text=label, showarrow=False, font=dict(size=12),
						xanchor="left", yanchor="middle")

	# Display the figure
	fig.show()
