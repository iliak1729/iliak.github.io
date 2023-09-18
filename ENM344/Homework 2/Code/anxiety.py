# Complete standalone Python script with data import and histograms based on age groups, degrees, and employment status
import csv
import matplotlib.pyplot as plt

# Initialize variables to store the data
age_satisfaction_dict = {}
degree_satisfaction_dict = {}
employment_satisfaction_dict = {}

# Read the CSV file using the csv library and ISO-8859-1 encoding
with open('Homework 2/anxiety.csv', 'r', encoding='ISO-8859-1') as csvfile:
    csvreader = csv.reader(csvfile)
    
    # Get the header to find the index of "SWL_T", "Age", "Degree", and "Employment"
    header = next(csvreader)
    swl_t_index = header.index("SWL_T")
    age_index = header.index("Age")
    degree_index = header.index("Degree")
    employment_index = header.index("Work")
    
    # Loop through the rows to collect data
    for row in csvreader:
        try:
            age = int(row[age_index])
            degree = row[degree_index]
            employment = row[employment_index]
            swl_t = float(row[swl_t_index])
            
            # Update the age_satisfaction_dict
            if age in age_satisfaction_dict:
                age_satisfaction_dict[age].append(swl_t)
            else:
                age_satisfaction_dict[age] = [swl_t]
                
            # Update the degree_satisfaction_dict
            if degree in degree_satisfaction_dict:
                degree_satisfaction_dict[degree].append(swl_t)
            else:
                degree_satisfaction_dict[degree] = [swl_t]
                
            # Update the employment_satisfaction_dict
            if employment in employment_satisfaction_dict:
                employment_satisfaction_dict[employment].append(swl_t)
            else:
                employment_satisfaction_dict[employment] = [swl_t]
                    
        except ValueError:
            # Skip rows where conversion to int or float fails
            continue

# Code for Age Groups Plot
sorted_ages = sorted(age_satisfaction_dict.keys())
age_group_intervals = list(range(min(sorted_ages), max(sorted_ages) + 1, 5))
age_group_satisfaction_dict = {}
for age in sorted_ages:
    for i in range(len(age_group_intervals) - 1):
        lower_bound = age_group_intervals[i]
        upper_bound = age_group_intervals[i + 1]
        if lower_bound <= age < upper_bound:
            age_group = f"{lower_bound}-{upper_bound - 1}"
            if age_group in age_group_satisfaction_dict:
                age_group_satisfaction_dict[age_group].extend(age_satisfaction_dict[age])
            else:
                age_group_satisfaction_dict[age_group] = age_satisfaction_dict[age]
            break
age_group_satisfaction_avg = {age_group: sum(swl_t_list) / len(swl_t_list) for age_group, swl_t_list in age_group_satisfaction_dict.items()}
sorted_age_groups = sorted(age_group_satisfaction_avg.keys(), key=lambda x: int(x.split('-')[0]))
sorted_avg_satisfaction_group = [age_group_satisfaction_avg[age_group] for age_group in sorted_age_groups]
plt.figure(figsize=(12, 6))
plt.bar(sorted_age_groups, sorted_avg_satisfaction_group, color='purple')
plt.xlabel('Age Group')
plt.ylabel('Average Satisfaction with Life')
plt.title('Average Satisfaction with Life as a Function of Age Groups')
plt.show()

# Code for Degree Plot
degree_satisfaction_avg = {degree: sum(swl_t_list) / len(swl_t_list) for degree, swl_t_list in degree_satisfaction_dict.items()}
sorted_degrees = sorted(degree_satisfaction_avg.keys())
sorted_avg_satisfaction_degree = [degree_satisfaction_avg[degree] for degree in sorted_degrees]
plt.figure(figsize=(15, 6))
plt.bar(sorted_degrees, sorted_avg_satisfaction_degree, color='green')
plt.xlabel('Degree')
plt.ylabel('Average Satisfaction with Life')
plt.title('Average Satisfaction with Life as a Function of Degree')
plt.xticks(rotation=45, ha='right')
plt.show()

# Code for Employment Status Plot
employment_satisfaction_avg = {employment: sum(swl_t_list) / len(swl_t_list) for employment, swl_t_list in employment_satisfaction_dict.items()}
sorted_employment = sorted(employment_satisfaction_avg.keys())
sorted_avg_satisfaction_employment = [employment_satisfaction_avg[employment] for employment in sorted_employment]
plt.figure(figsize=(8, 6))
plt.bar(sorted_employment, sorted_avg_satisfaction_employment, color='orange')
plt.xlabel('Employment Status')
plt.ylabel('Average Satisfaction with Life')
plt.title('Average Satisfaction with Life as a Function of Employment Status')
plt.xticks(rotation=45, ha='right')
plt.show()

# Run the complete standalone Python script with data import and histograms based on age groups, degrees, work status, and age groups with separate columns for each level of education

# Initialize variables to store the data
age_degree_satisfaction_dict = {}

# Read the CSV file using the csv library and ISO-8859-1 encoding
with open('Homework 2/anxiety.csv', 'r', encoding='ISO-8859-1') as csvfile:
    csvreader = csv.reader(csvfile)
    
    # Get the header to find the index of "SWL_T", "Age", and "Degree"
    header = next(csvreader)
    swl_t_index = header.index("SWL_T")
    age_index = header.index("Age")
    degree_index = header.index("Degree")
    
    # Loop through the rows to collect data
    for row in csvreader:
        try:
            age = int(row[age_index])
            degree = row[degree_index]
            swl_t = float(row[swl_t_index])
            
            # Update the age_degree_satisfaction_dict
            age_degree_key = (age, degree)
            if age_degree_key in age_degree_satisfaction_dict:
                age_degree_satisfaction_dict[age_degree_key].append(swl_t)
            else:
                age_degree_satisfaction_dict[age_degree_key] = [swl_t]
                    
        except ValueError:
            # Skip rows where conversion to int or float fails
            continue

# Code for Age Groups Plot with Separate Columns for Each Level of Education
sorted_ages = sorted(list(set([key[0] for key in age_degree_satisfaction_dict.keys()])))
age_group_intervals = list(range(min(sorted_ages), max(sorted_ages) + 1, 5))
age_degree_group_satisfaction_dict = {}
for age_degree_key, swl_t_list in age_degree_satisfaction_dict.items():
    age, degree = age_degree_key
    for i in range(len(age_group_intervals) - 1):
        lower_bound = age_group_intervals[i]
        upper_bound = age_group_intervals[i + 1]
        if lower_bound <= age < upper_bound:
            age_group = f"{lower_bound}-{upper_bound - 1}"
            age_degree_group_key = (age_group, degree)
            if age_degree_group_key in age_degree_group_satisfaction_dict:
                age_degree_group_satisfaction_dict[age_degree_group_key].extend(swl_t_list)
            else:
                age_degree_group_satisfaction_dict[age_degree_group_key] = swl_t_list
            break

age_degree_group_satisfaction_avg = {age_degree_group_key: sum(swl_t_list) / len(swl_t_list) for age_degree_group_key, swl_t_list in age_degree_group_satisfaction_dict.items()}
sorted_age_degree_groups = sorted(list(set([key[0] for key in age_degree_group_satisfaction_avg.keys()])), key=lambda x: int(x.split('-')[0]))
sorted_degrees = sorted(list(set([key[1] for key in age_degree_group_satisfaction_avg.keys()])))

plt.figure(figsize=(15, 8))

for degree in sorted_degrees:
    sorted_avg_satisfaction_age_degree_group = [age_degree_group_satisfaction_avg.get((age_group, degree), 0) for age_group in sorted_age_degree_groups]
    plt.bar(sorted_age_degree_groups, sorted_avg_satisfaction_age_degree_group, alpha=0.7, label=degree)

plt.xlabel('Age Group')
plt.ylabel('Average Satisfaction with Life')
plt.title('Average Satisfaction with Life as a Function of Age Groups (Separated by Degree)')
plt.legend(title='Degree')
plt.show()

# Initialize a variable to store the data for the stacked bar graph with 1-year age groups
age_degree_grouped_stacked_1yr_data = {}

# Define degree groups as they are reported in the data file, accounting for variations in naming
degree_groups_for_plotting = ['Bachelor', 'High School Diploma', 'Masters', 'None', 'Ph.D.']

# Function to categorize degrees
def categorize_degree(degree):
    if 'Bachelor' in degree:
        return 'Bachelor'
    elif 'High School' in degree:
        return 'High School Diploma'
    elif 'Master' in degree:
        return 'Masters'
    elif 'None' in degree:
        return 'None'
    elif 'Ph.D.' in degree:
        return 'Ph.D.'
    else:
        return 'Other'

# Organize the data for the stacked bar graph with 1-year age groups
for age_degree_key, swl_t_list in age_degree_satisfaction_dict.items():
    age, degree = age_degree_key
    age_group_1yr = str(age)
    
    # Use function to categorize degrees
    degree_group = categorize_degree(degree)
    
    age_degree_grouped_stacked_1yr_key = (age_group_1yr, degree_group)
    if age_degree_grouped_stacked_1yr_key in age_degree_grouped_stacked_1yr_data:
        age_degree_grouped_stacked_1yr_data[age_degree_grouped_stacked_1yr_key] += len(swl_t_list)
    else:
        age_degree_grouped_stacked_1yr_data[age_degree_grouped_stacked_1yr_key] = len(swl_t_list)

# Prepare the data for plotting with 1-year age groups
age_groups_1yr_for_plotting = sorted(list(set([key[0] for key in age_degree_grouped_stacked_1yr_data.keys()])))

# Initialize the plot
plt.figure(figsize=(20, 8))
bottoms_1yr = [0] * len(age_groups_1yr_for_plotting)

# Define contrasting colors
contrasting_colors = ['b', 'g', 'r', 'c', 'm', 'y']

# Loop through each degree and add a layer to the stacked bar graph
for idx, degree in enumerate(degree_groups_for_plotting):
    color = contrasting_colors[idx % len(contrasting_colors)]
    counts_for_this_degree_1yr = [age_degree_grouped_stacked_1yr_data.get((age_group, degree), 0) for age_group in age_groups_1yr_for_plotting]
    plt.bar(age_groups_1yr_for_plotting, counts_for_this_degree_1yr, color=color, label=degree, bottom=bottoms_1yr)
    bottoms_1yr = [bottom + count for bottom, count in zip(bottoms_1yr, counts_for_this_degree_1yr)]

# Add labels, title, and legend
plt.xlabel('Age (1-year intervals)')
plt.ylabel('Count of Individuals')
plt.title('Stacked Bar Graph of Individual Ages vs Degree Level')
plt.legend(title='Degree Level')
plt.xticks(rotation=90)
plt.show()