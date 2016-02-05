function [ U ] = stage1( A )
%stage1 Converts a real, square matrix to upper echelon form
%   Converts a real, square matrix to upper echelon form, if possible

tic

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

% Put largest cell values onto diagonal (we'll use the diagonal for our 
% pivot later)
for column_num = 1:column_count-1
    row_inspecting = column_num;
    row_with_pivot = column_num;
    
    while row_inspecting <= row_count
        if U(row_inspecting, column_num) > U(row_with_pivot, column_num)
            row_with_pivot = row_inspecting;
        end
        
        row_inspecting = row_inspecting + 1;
    end
    
    % Swap the row with the largest pivot for the current column into place
    if(column_num == row_with_pivot)
       % Do nothing... It's already in the right place 
    else 
        temp = U(row_with_pivot, :);
        U(row_with_pivot, :) = U(column_num, :);
        U(column_num, :) = temp;
    end
end

% For each column, 1->n
for column_inspecting = 1:column_count
    % For each row from 2->n
    for row_inspecting = 2:row_count
        
        % Make mutator/pivot be U(n,n) (because that's the largest pivot 
        % as we set it to be so earlier 
        mutator_row_number = column_inspecting;
        
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
full_rank = 1;
rank_check_row = row_count;

while rank_check_row >= 1;
    num_zeroes_in_row = 0;
    
    for rank_check_column = 1:column_count
        if U(rank_check_row, rank_check_column) == 0
            num_zeroes_in_row = num_zeroes_in_row + 1; 
        end
    end
    
    if num_zeroes_in_row == column_count
       fprintf('Matrix does not have full rank \n')
       
       full_rank = 0;
       
       break;
    end
    
    rank_check_row = rank_check_row - 1;
end

if full_rank == 1
    fprintf('Matrix has full rank.\n')
end

toc

end

