const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};
export default Link;
