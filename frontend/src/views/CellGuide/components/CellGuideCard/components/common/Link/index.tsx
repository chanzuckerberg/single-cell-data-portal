import { ReactElement } from "react";
import { StyledLink } from "./style";

const Link = ({
  label,
  url,
  onClick,
}: {
  label: string | ReactElement;
  url: string;
  onClick?: () => void;
}) => {
  return (
    <StyledLink href={url} target="_blank" onClick={onClick}>
      {label}
    </StyledLink>
  );
};
export default Link;
