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

    % For each column, 1->n
    for column_inspecting = 1:column_count
        % Make mutator/pivot be U(n,n) (because that's the largest pivot 
        % as we set it to be so earlier 
        mutator_row_number = column_inspecting;

        % Pick the largest pivot/mutator at or below the existing 
        % mutator_row_number
        for row_inspecting = mutator_row_number:row_count
            if U(row_inspecting, column_inspecting) > U(mutator_row_number, column_inspecting)
               mutator_row_number = row_inspecting;
            end
        end

        % Swap the mutator row so the pivot is on the diagonal
        if mutator_row_number == column_inspecting
            % Do nothing... It's already on the diagonal
        else 
           temp = U(mutator_row_number, :);
           U(mutator_row_number, :) = U(column_inspecting, :);
           U(column_inspecting, :) = temp;
        end

        % Set the mutator_row_number to be the diagonal again
        mutator_row_number = column_inspecting;

        % For each row from 2->n
        for row_inspecting = 2:row_count
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

                % Find out what we have to divide our mutator row by in order
                % to create a 0 in the row/column inspecting
                multiplication_factor = cell_value_in_current_row / cell_value_in_mutator_row;

                % Create the transformed mutator row, using our multiplication
                % factor
                mutator_row = U(mutator_row_number, :) .* multiplication_factor;

                % Subtract our mutator row from the row we're inspecting
                U(row_inspecting, :) = U(row_inspecting, :) - mutator_row;
            end
        end
    end


    % Check that no rows contain only 0's (if they do, it's not full rank)
    full_rank = 1;
    rank_check_row = row_count;

    % Check the last row first, because it's the one most likely contain
    % only 0's
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
end

