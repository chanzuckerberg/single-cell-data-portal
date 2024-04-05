import { ReactElement } from "react";
import { StyledLink } from "./style";

const Link = ({
  label,
  url,
  onClick,
  dataTestId,
}: {
  label: string | ReactElement;
  url: string;
  onClick?: () => void;
  dataTestId?: string;
}) => {
  return (
    <StyledLink
      data-testid={dataTestId}
      href={url}
      target="_blank"
      onClick={onClick}
    >
      {label}
    </StyledLink>
  );
};
export default Link;
