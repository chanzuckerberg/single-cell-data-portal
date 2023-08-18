import { ReactElement } from "react";
import { StyledLink } from "./style";

const Link = ({
  label,
  url,
}: {
  label: string | ReactElement;
  url: string;
}) => {
  return (
    <StyledLink href={url} target="_blank">
      {label}
    </StyledLink>
  );
};
export default Link;
