function [ U ] = stage1( A )
%stage1 Converts a real, square matrix to upper echelon form
%   Converts a real, square matrix to upper echelon form, if possible


% Check that it's square
dimensions = size(A);
if dimensions(1,1) ~= dimensions(1,2)
    error('Input matrix is not square')
end

% Check that it's real
for index = 1:numel(A)
    if A(index)*A(index) < 0
        error('Input matrix has imaginary number(s)')
    end
end

row_count = dimensions(1,1);
column_count = dimensions(1,2);
U = A;

% Move rows with 0s to the end, such that the higher the column number of 
% the cell conaining a 0, the further down the row is in the matrix
for column_num = 1:column_count
    rows_inspected = 0;
    row_num = 1;
    while row_num < row_count
        % If the current cell being inspected is 0
        if U(row_num, column_num) == 0
            % Save the row to a temp vector
            temp = U(row_num, :);
            
            % Shift all rows below it up one
            for row_below = row_num + 1:row_count
                U(row_below - 1, :) = U(row_below, :);
            end
            
            % Shift the current row to the bottom
            U(row_count, :) = temp;
                        
            % If there are lots of zeroes, we may end up continously moving
            % the current row to the bottom in an infinite loop.  We can
            % avoid this with the following if/break
            if rows_inspected == row_count - 1
               break; 
            end      
        else          
            row_num = row_num + 1;
        end
        
        rows_inspected = rows_inspected + 1;
    end
end

% For each column, 1->n
for column_inspecting = 1:column_count
    % For each row from 2->n
    for row_inspecting = 2:row_count
        
        % Make mutator/pivot row be the biggest above the current row being inspected 
        mutator_row_number = column_inspecting;
        for mutator_row_number_search = row_inspecting - 1:1
            if abs(U(mutator_row_number_search, column_inspecting)) > abs(U(mutator_row_number, column_inspecting))
                mutator_row_number = mutator_row_number_search;
            end
        end
        
        % If the column we're looking at is under the diagonal
        if column_inspecting < row_inspecting
            % Make that cell 0 using the mutator row.  As long as we use
            % a row above the current row, we'll never unset a 0
            % in a previous column
            cell_value_in_current_row = U(row_inspecting, column_inspecting);
            
            % If it's already 0, we don't have to do anything
            if cell_value_in_current_row == 0
                continue
            end
            
            cell_value_in_mutator_row = U(mutator_row_number, column_inspecting);
            
            % We can't divide by 0... This state *should* never happen
            if cell_value_in_mutator_row == 0
                error('Not enough info to convert to upper echelon form.')
            end
            
            multiplication_factor = cell_value_in_current_row / cell_value_in_mutator_row;
            mutator_row = U(mutator_row_number, :) .* multiplication_factor;
            
            U(row_inspecting, :) = U(row_inspecting, :) - mutator_row;
        end
    end
end


% Check that no rows contain only 0's (if they do, it's not full rank)
if U(row_count, column_count) == 0
    fprintf('Matrix does not have full rank.')
else
    fprintf('Matrix has full rank.')
end

end

