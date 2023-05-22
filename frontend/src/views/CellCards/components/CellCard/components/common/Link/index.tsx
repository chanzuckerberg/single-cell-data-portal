import { StyledLink } from "./style";

const Link = ({ label, url }: { label: string; url: string }) => {
  return (
    <StyledLink href={url} target="_blank">
      {label}
    </StyledLink>
  );
};
export default Link;
