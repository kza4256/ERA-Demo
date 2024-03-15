# ERA-Demo
This MATLAB script uses ERA to estimate a dynamical system's state-space model from input-output data. It includes functions to compute the Hankel matrix, perform SVD, and truncate SVD for model reduction. The inferred model is then used to predict the system's output, and the Frobenius error between the predicted and true outputs is calculated.

# To run the ERA script, follow these steps:

   1. Clone the Repository: Clone the repository containing the code to your local machine.
   2. Open MATLAB: Open MATLAB on your machine.
   3. Navigate to the Code Directory: Use the cd command to navigate to the directory where you cloned the repository.
   4. Run the Script: Run the script by typing ErrAnalysis(interp_points, m, ns) in the MATLAB command window. Replace interp_points, m, and ns with the desired values for interpolation points, m, and ns.
   5. View Results: The script will compute the ERA and visualize the true and predicted outputs in a MATLAB figure window.
   6. Adjust Parameters: You can adjust the parameters and input data in the script to analyze different dynamical systems and datasets.
   7. Additional Information: For more information about the ERA algorithm and the functions used in the script, refer to the comments within the script and the function descriptions in the README.

# Note:
Ensure that you have the required data files (ch01PreCO2_cps_values, ch01PreCO2_cps_times, etc.) in the same directory as the script or provide the correct path to these files in the script.
