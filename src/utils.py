
def list_split(input_list, num_splits):
    """ Split a list into multiple sub-lists. """

    if num_splits > len(input_list):
        raise ValueError("Cannot split a list with more splits than its actual size.")

    # calculate the approximate size of each sublist
    avg_size = len(input_list) // num_splits
    remainder = len(input_list) % num_splits

    # initialize variables
    start = 0
    end = avg_size
    sublists = []

    for i in range(num_splits):
        # adjust sublist size for the remainder
        if i < remainder:
            end += 1

        # create a sublist and add it to the result
        sublist = input_list[start:end]
        sublists.append(sublist)

        # update the start and end indices for the next sublist
        start = end
        end += avg_size

    return sublists


def divide_grid(grid, n_total, task_id):
    """ Divide a grid into n_total tasks and return the task corresponding to task_id.

    :param grid: list of arguments for each call
    :param n_total: number of workers
    :param task_id: id of this worker
    """
    # total number of tasks exceeds the grid size
    if n_total > len(grid):
        n_total = len(grid)
    # task id invalid: trying to complete null tasks
    if task_id > n_total:
        return None
    # task_id corresponds to a chunk of the grid
    return list_split(grid, n_total)[task_id]
