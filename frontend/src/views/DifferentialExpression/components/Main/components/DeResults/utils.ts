const parseExpressions = (expression: string) => {
  const expressions = expression
    .split(",")
    .map((expr) => {
      const match = expr.match(/([<>]=?)\|?(-?\d*\.?\d+(?:e[+-]?\d+)?)\|?/);
      if (!match) return null;

      const [, operator, valueStr] = match;
      const value = parseFloat(valueStr);
      const isAbsolute = expr.includes("|");
      return { operator, value, isAbsolute };
    })
    .filter((expr) => expr !== null);

  return expressions.length > 0 ? expressions : null;
};

export const applyFilter = (resultValue: number, filterExpression: string) => {
  const parsedExpressions = parseExpressions(filterExpression);
  if (!parsedExpressions) return true;
  const filteredParsedExpressions = parsedExpressions.filter(
    (expr) => expr !== null
  ) as {
    operator: string;
    value: number;
    isAbsolute: boolean;
  }[];
  return filteredParsedExpressions.every(({ operator, value, isAbsolute }) => {
    const targetValue = isAbsolute ? Math.abs(resultValue) : resultValue;
    switch (operator) {
      case "<":
        return targetValue < value;
      case "<=":
        return targetValue <= value;
      case ">":
        return targetValue > value;
      case ">=":
        return targetValue >= value;
      default:
        return true;
    }
  });
};
