import re
import subprocess


def run_mypy():
    result = subprocess.run(["mypy", "--ignore-missing-imports", "."], capture_output=True, text=True)
    return result.stdout


def annotate_errors(output):
    error_pattern = re.compile(r"([^:]+:\d+): error:.*$")
    for line in output.splitlines():
        match = error_pattern.match(line)
        if match:
            file_and_line = match.group(1)
            file_path, line_number = file_and_line.split(":")
            if file_path.endswith(".py"):
                with open(file_path, "r+") as file:
                    lines = file.readlines()
                    line_number = int(line_number)
                    if line_number <= len(lines):
                        line = lines[line_number - 1].rstrip()  # Remove trailing newline
                        if "# type: ignore" not in line:
                            comment_index = line.find("#")
                            if comment_index >= 0:
                                line_parts = [line[:comment_index], line[comment_index:]]
                            else:
                                line_parts = [line, ""]
                            updated_line = f"{line_parts[0]}  # type: ignore{line_parts[1]}\n"
                            lines[line_number - 1] = updated_line
                            file.seek(0)
                            file.writelines(lines)


if __name__ == "__main__":
    mypy_output = run_mypy()
    annotate_errors(mypy_output)
