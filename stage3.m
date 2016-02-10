function [ x ] = stage3( A, b )
%stage1 Solves for x in equation Ax = b
%   Returns a solution vector given a coefficient matrix and
%   right-hand-side values
%
%   By Jonathan Sterling - u1276062

    % Check that A is square
    dimensions = size(A);
    if dimensions(1,1) ~= dimensions(1,2)
        error('Input coefficient matrix is not square')
    end
    
    % Check that A is not empty
    if(dimensions(1,1) == 0 || dimensions(1,2) == 0)
        error('Input coefficient matrix is empty')
    end

    % Check that b is the right size (also checks that it's not empty,
    % because A is not empty
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

    % Get our the row and column counts of A, then merge A and b into an
    % augmented matrix called AUG
    row_count = dimensions(1,1);
    column_count = dimensions(1,2);
    AUG = [A,b];

    % For each column, 1->n
    for column_inspecting = 1:column_count
        % Initialize mutator row number to the top possible column...
        mutator_row_number = column_inspecting;

        % Pick the largest pivot/mutator at or below the existing 
        % mutator_row_number
        for row_inspecting = column_inspecting:row_count
            if abs(AUG(row_inspecting, column_inspecting)) > abs(AUG(mutator_row_number, column_inspecting))
               mutator_row_number = row_inspecting;
            end
        end

        % Swap the mutator row so the pivot is on the diagonal
        if mutator_row_number == column_inspecting
            % Do nothing... It's already on the diagonal
        else 
           % Put the mutator row into a temporary vector
           temp_A = AUG(mutator_row_number, :);

           % Move row on the diagonal to the mutator row
           AUG(mutator_row_number, :) = AUG(column_inspecting, :);

           % Move the mutator row (stored in temp_A) so that the pivot is
           % on the diagonal
           AUG(column_inspecting, :) = temp_A;
        end

        % Set the mutator_row_number to be the diagonal again
        mutator_row_number = column_inspecting;

        % For each row from under the mutator row
        for row_inspecting = (mutator_row_number + 1):row_count     
            % Make that cell 0 using the mutator row.
            cell_value_in_current_row = AUG(row_inspecting, column_inspecting);

            % If it's already 0, we don't have to do anything
            if cell_value_in_current_row == 0
                continue
            end

            cell_value_in_mutator_row = AUG(mutator_row_number, column_inspecting);

            % We can't divide by 0... This state *should* never happen
            if cell_value_in_mutator_row == 0
                error('Not enough info to convert to upper echelon form.')
            end

            % Find out what we have to divide our mutator row by in order
            % to create a 0 in the row/column inspecting
            multiplication_factor = cell_value_in_current_row / cell_value_in_mutator_row;

            % Create the transformed mutator rows, using our multiplication
            % factor
            mutator_row_A = AUG(mutator_row_number, :) .* multiplication_factor;

            % Subtract our mutator row from the row we're inspecting
            AUG(row_inspecting, :) = AUG(row_inspecting, :) - mutator_row_A;
        end
    end

    % If the last cell is 0, the last row is 0, and the matrix is rank
    % deficient
    rank_deficient = 0;
    if AUG(row_count, column_count) == 0 
        rank_deficient = 1;
    end
    
    if rank_deficient == 1
        error('Cannot solve.  Matrix does not have full rank')
    else
        fprintf('Matrix has full rank. Solving via back substitution.\n')
    end

    % Solve from last row back to first...
    x = zeros(row_count,1);
    row_solving = row_count;

    while row_solving >= 1
        % If we're solving for X1 in (#1*X1 + #2*X2 + #3*X3 = 4),
        % calculate #2*X2 and #3*X3, subtract that value from both sides of
        % the equation, then divide both sides by #1 to get the value of X1
        
        value_of_row_after_unknown = 0;
        column_index = row_solving + 1;
        
        while column_index <= column_count
            value_of_row_after_unknown = value_of_row_after_unknown + (AUG(row_solving, column_index) * x(column_index));

            column_index = column_index + 1;
        end

        x(row_solving) = (AUG(row_solving, column_count + 1) - value_of_row_after_unknown) / AUG(row_solving, row_solving);

        row_solving = row_solving - 1;
    end
end