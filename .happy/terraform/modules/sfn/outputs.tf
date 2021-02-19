output step_function_arn {
  value       = aws_sfn_state_machine.state_machine.id
  description = "ARN for the step function definition"
}
