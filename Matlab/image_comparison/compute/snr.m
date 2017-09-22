function [snrvalue, snrvalue_sigma] = snr(image_data, image_solution, sigma_value)

    [height,width] = size(image_data);

    vector_data = reshape(image_data, height*width, 1);
    vector_solution = reshape(image_solution, height*width, 1);

    std_image_solution = std(vector_solution);
    std_image_data = std(vector_data);

    snrvalue = std_image_solution/std_image_data;
    snrvalue_sigma = std_image_solution/sigma_value;

end