function [ x ] = stage3( A, b )
%stage1 Solves for x in equation Ax = b
%   Returns a solution vector given a coefficient matrix and
%   right-hand-side values.  Uses back substitution, partial pivoting, and
%   scaling

tic

% Check that it's square
dimensions = size(A);
if dimensions(1,1) ~= dimensions(1,2)
    error('Input coefficient matrix is not square')
end

% Check that b is the right size
if length(b) ~= dimensions(1,1)
   error('Input values vector is not the correct size') 
end

% Check that A contains only real numbers
for index = 1:numel(A)
    if A(index)*A(index) < 0
        error('Input coefficient matrix has imaginary number(s)')
    end
end

% Check that b contains only real numbers
for index = 1:length(b)
   if b(index)*b(index) < 0
       error('Input values vector has imaginary number(s)');
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
        temp_A = U(row_with_pivot, :);
        temp_b = b(row_with_pivot);
        
        U(row_with_pivot, :) = U(column_num, :);
        b(row_with_pivot) = b(column_num);
        
        U(column_num, :) = temp_A;
        b(column_num) = temp_b;
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
            mutator_row_A = U(mutator_row_number, :) .* multiplication_factor;
            mutator_row_b = b(mutator_row_number) .* multiplication_factor;
            
            U(row_inspecting, :) = U(row_inspecting, :) - mutator_row_A;
            b(row_inspecting) = b(row_inspecting) - mutator_row_b;
        end
    end
end


% Check that no rows contain only 0's (if they do, it's not full rank)
rank_check_row = row_count;

while rank_check_row >= 1
    num_zeroes_in_row = 0;
    
    for rank_check_column = 1:column_count
        if U(rank_check_row, rank_check_column) == 0
            num_zeroes_in_row = num_zeroes_in_row + 1; 
        end
    end
    
    if num_zeroes_in_row == column_count
       error('Cannot solve.  Matrix does not have full rank')
    end
    
    rank_check_row = rank_check_row - 1;
end

fprintf('Matrix has full rank. Solving via back substitution.\n')

% Solve from last row back to first...
x = zeros(row_count,1);
row_solving = row_count;

while row_solving >= 1
    value_of_row_after_unknown = 0;
    
    column_index = row_solving;
    while column_index <= column_count
        value_of_row_after_unknown = value_of_row_after_unknown + (U(row_solving, column_index) * x(column_index));
        
        column_index = column_index + 1;
    end
    
    x(row_solving) = (b(row_solving) - value_of_row_after_unknown) / U(row_solving, row_solving);
    
    row_solving = row_solving - 1;
end

toc

end

