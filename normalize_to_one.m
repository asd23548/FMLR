function output_array = normalize_to_one(input_array,col_index)
% Function to normalize an array
temp_array = input_array(:,col_index);
sum_array = sum(temp_array,2);
output_array = temp_array./sum_array;
