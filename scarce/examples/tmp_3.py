import numpy as np
import matplotlib.pyplot as plt

t_steps = 1000
t = np.linspace(0, 3 * np.pi, t_steps)
y = np.array([-np.log(t / 10. + 0.1),
              t,
              np.log(t / 10. + 0.1) + t**2]).T

y -= y.min(axis=0)
y += 0.1

y /= y.max(axis=0)

q_max = 1.

# Result arrays
n_steps = 100
t_step = np.zeros((n_steps, y.shape[1]))
y_step = np.zeros_like(t_step)
i_step = np.zeros(y.shape[1], dtype=np.int)  # Result arrays storage indeces
next_step = np.zeros_like(i_step)

step_size = t_steps / n_steps

# print 't.shape', t.shape
# print 'y.shape', y.shape
# print 't_step.shape', t_step.shape

old_val = None
old_x = None

max_step_size = t_steps / n_steps * 8


def cal_step_size(t, y, dt, dy, q_max, i_step, n_steps):
    ''' Calculates the step size from the actual data

    It is assumed that the actual slope stays constant.
    The remaining time distance is calculated and the step size is
    adjusted to fit the remaining steps.
    '''

    step_size = np.zeros_like(dy)

    # All steps - 1 are used, mark as done
    sel_done = n_steps - i_step - 1 <= 0

    # All storage spaces used up
    if step_size[~sel_done].size == 0:
        return step_size.astype(np.int)

    # Set next step to NaN is storing is done
    step_size[sel_done] = np.NaN

    # Calculate slope for linear data extrapolation
    dy_dt = dy[~sel_done] / dt

    # Calculate remaining x distance to cover

    # Case: increasing function
    t_max_exp = (q_max - y[~sel_done]) / dy_dt + t
    # Case: increasing function
    t_max_exp[dy_dt < 0] = ((0 - y[~sel_done]) / dy_dt + t)[dy_dt < 0]
    # Remaining time to cover
    t_left = t_max_exp - t

    # Calculate the step size
    step_size[~sel_done] = t_left / dt / (n_steps - i_step[~sel_done] - 1)

    # Limit step size to 5 times the distance of a constant
    # step size. Needed for slope direction changing funcctions
    sel_limit = step_size[~sel_done] > max_step_size
    step_size[np.logical_and(~sel_done, sel_limit)] = max_step_size

    step_size = step_size.astype(np.int)
    return step_size


def store_if_needed(i, next_step, t, y, dt, dy, q_max, i_step, n_steps):
    ''' Checks if charge carriers value needs to be stored and returns
        next storage time step
    '''
    # Special case at the beginning: calculate 2nd storage step
    # from first 2 time points
    if i == 1:
        next_step += cal_step_size(t[i], y[i], dt, dy, q_max, i_step, n_steps)
        return

    # Select charges that need storing
    sel = i == next_step

    # No carrier needs storing
    if i_step[sel].size == 0:
        return

    # Set data of actual time step
    t_step[i_step[sel], sel] = t[i]
    y_step[i_step[sel], sel] = y[i, sel]

    # Increase storage hists indeces
    i_step[sel] += 1

    # Set new storage index
    if dy is not None:  # Happens on first storage at t = 0
        next_step[sel] += cal_step_size(t[i], y[i, sel], dt, dy[sel], q_max, i_step[sel], n_steps)
        print 'new next_step', next_step

dt = np.diff(t)[0]
dy = None

for i, val in enumerate(y):
    if old_val is not None:
        dy = val - old_val

    store_if_needed(i, next_step, t, y, dt, dy, q_max, i_step, n_steps)

    old_val = val
    old_x = t[i]

# Force last save point to last data point
t_step[-1] = t[-1]
y_step[-1] = y[-1]


plt.plot(t, y)
plt.plot(t_step, y_step, 'o')
plt.show()
